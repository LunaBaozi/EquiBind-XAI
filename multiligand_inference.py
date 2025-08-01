#!/usr/bin/env python
import argparse
import sys
import dgl

from copy import deepcopy

import os
import numpy as np

from rdkit import Chem
from rdkit.Geometry import Point3D

from commons.geometry_utils import rigid_transform_Kabsch_3D, get_torsions, get_dihedral_vonMises, apply_changes
from commons.process_mols import get_rec_graph, get_receptor_inference

#from train import load_model

from commons.utils import seed_all

import yaml

from datasets.custom_collate import *  # do not remove
from models import *  # do not remove
from torch.nn import *  # do not remove
from torch.optim import *  # do not remove
from commons.losses import *  # do not remove
from torch.optim.lr_scheduler import *  # do not remove
from torch.utils.data import DataLoader


# turn on for debugging C code like Segmentation Faults
import faulthandler
from datasets import multiple_ligands

faulthandler.enable()

from models.equibind import EquiBind

def parse_arguments(arglist = None):
    p = argparse.ArgumentParser()    
    p.add_argument("-l", "--ligands_sdf", type=str, help = "A single sdf file containing all ligands to be screened when running in screening mode")
    p.add_argument("-r", "--rec_pdb", type = str, help = "The receptor to dock the ligands in --ligands_sdf against")
    p.add_argument('-o', '--output_directory', type=str, default=None, help='path where to put the predicted results')
    p.add_argument('--config', type=argparse.FileType(mode='r'), default=None)
    p.add_argument('--checkpoint', '--model', dest = "checkpoint",
                   type=str, help='path to .pt file containing the model used for inference. '
                   'Defaults to runs/flexible_self_docking/best_checkpoint.pt in the same directory as the file being run')
    p.add_argument('--train_args', type = str, help = "Path to a yaml file containing the parameters that were used to train the model. "
                    "If not supplied, it is assumed that a file named 'train_arguments.yaml' is located in the same directory as the model checkpoint")
    p.add_argument('--no_skip', dest = "skip_in_output", action = "store_false", help = 'skip input files that already have corresponding folders in the output directory. Used to resume a large interrupted computation')
    p.add_argument('--batch_size', type=int, default=8, help='samples that will be processed in parallel')
    p.add_argument("--n_workers_data_load", type = int, default = 0, help = "The number of cores used for loading the ligands and generating the graphs used as input to the model. 0 means run in correct process.")
    p.add_argument('--use_rdkit_coords', action="store_true", help='override the rkdit usage behavior of the used model')
    p.add_argument('--device', type=str, default='cuda', help='What device to train on: cuda or cpu')
    p.add_argument('--seed', type=int, default=1, help='seed for reproducibility')
    p.add_argument('--num_confs', type=int, default=1, help='num_confs if using rdkit conformers')
    p.add_argument("--lig_slice", help = "Run only a slice of the provided ligand file. Like in python, this slice is HALF-OPEN. Should be provided in the format --lig_slice start,end")
    p.add_argument("--lazy_dataload", dest = "lazy_dataload", action="store_true", default = None, help = "Turns on lazy dataloading. If on, will postpone rdkit parsing of each ligand until it is requested.")
    p.add_argument("--no_lazy_dataload", dest = "lazy_dataload", action="store_false", default = None, help = "Turns off lazy dataloading. If on, will postpone rdkit parsing of each ligand until it is requested.")
    p.add_argument("--no_run_corrections", dest = "run_corrections", action = "store_false", help = "possibility of turning off running fast point cloud ligand fitting")

    cmdline_parser = deepcopy(p)
    args = p.parse_args(arglist)
    clear_defaults = {key: argparse.SUPPRESS for key in args.__dict__}
    cmdline_parser.set_defaults(**clear_defaults)
    cmdline_parser._defaults = {}
    cmdline_args = cmdline_parser.parse_args(arglist)

    return p.parse_args(arglist), set(cmdline_args.__dict__.keys())

def get_default_args(args, cmdline_args):
    if args.config:
        config_dict = yaml.load(args.config, Loader=yaml.FullLoader)
        arg_dict = args.__dict__
        for key, value in config_dict.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].append(v)
            else:
                arg_dict[key] = value
        args.config = args.config.name
    else:
        config_dict = {}
    
    if args.checkpoint is None:
        args.checkpoint = os.path.join(os.path.dirname(__file__), "runs/flexible_self_docking/best_checkpoint.pt")
    
    config_dict['checkpoint'] = args.checkpoint
    # overwrite args with args from checkpoint except for the args that were contained in the config file or provided directly in the commandline
    arg_dict = args.__dict__

    if args.train_args is None:
        with open(os.path.join(os.path.dirname(args.checkpoint), 'train_arguments.yaml'), 'r') as arg_file:
            checkpoint_dict = yaml.load(arg_file, Loader=yaml.FullLoader)
    else:
        with open(args.train_args, 'r') as arg_file:
            checkpoint_dict = yaml.load(arg_file, Loader=yaml.FullLoader)

    for key, value in checkpoint_dict.items():
        if (key not in config_dict.keys()) and (key not in cmdline_args):
            if isinstance(value, list) and key in arg_dict:
                for v in value:
                    arg_dict[key].append(v)
            else:
                arg_dict[key] = value
    args.model_parameters['noise_initial'] = 0
    return args

def load_rec_and_model(args):
    device = torch.device("cuda:0" if torch.cuda.is_available() and args.device == 'cuda' else "cpu")
    print(f"device = {device}")
    # sys.exit()
    checkpoint = torch.load(args.checkpoint, map_location=device)
    dp = args.dataset_params

    model = EquiBind(device = device, lig_input_edge_feats_dim = 15, rec_input_edge_feats_dim = 27, **args.model_parameters)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.to(device)
    model.eval()

    rec_path = args.rec_pdb
    rec, rec_coords, c_alpha_coords, n_coords, c_coords = get_receptor_inference(rec_path)
    rec_graph = get_rec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords,
                                use_rec_atoms=dp['use_rec_atoms'], rec_radius=dp['rec_graph_radius'],
                                surface_max_neighbors=dp['surface_max_neighbors'],
                                surface_graph_cutoff=dp['surface_graph_cutoff'],
                                surface_mesh_cutoff=dp['surface_mesh_cutoff'],
                                c_alpha_max_neighbors=dp['c_alpha_max_neighbors'])

    return rec_graph, model

# def run_batch(model, ligs, lig_coords, lig_graphs, rec_graphs, geometry_graphs, true_indices):
#     try:
#         predictions = model(lig_graphs, rec_graphs, geometry_graphs)[0]
#         out_ligs = ligs
#         out_lig_coords = lig_coords
#         names = [lig.GetProp("_Name") for lig in ligs]
#         successes = list(zip(true_indices, names))
#         failures = []
#     except AssertionError:
#         lig_graphs, rec_graphs, geometry_graphs = (dgl.unbatch(lig_graphs),
#         dgl.unbatch(rec_graphs), dgl.unbatch(geometry_graphs))
#         predictions = []
#         out_ligs = []
#         out_lig_coords = []
#         successes = []
#         failures = []
#         for lig, lig_coord, lig_graph, rec_graph, geometry_graph, true_index in zip(ligs, lig_coords, lig_graphs, rec_graphs, geometry_graphs, true_indices):
#             try:
#                 output = model(lig_graph, rec_graph, geometry_graph)
#             except AssertionError as e:
#                 failures.append((true_index, lig.GetProp("_Name")))
#                 print(f"Failed for {lig.GetProp('_Name')}")
#             else:
#                 out_ligs.append(lig)
#                 out_lig_coords.append(lig_coord)
#                 predictions.append(output[0][0])
#                 successes.append((true_index, lig.GetProp("_Name")))
#     assert len(predictions) == len(out_ligs)
#     return out_ligs, out_lig_coords, predictions, successes, failures

# Support for confidence score extraction
def run_batch(model, ligs, lig_coords, lig_graphs, rec_graphs, geometry_graphs, true_indices):
    try:
        # Get model output
        model_output = model(lig_graphs, rec_graphs, geometry_graphs)
        
        # Handle different model output formats
        if isinstance(model_output, tuple) and len(model_output) >= 6:
            predictions, ligs_keypts, recs_keypts, rotations, translations, geom_losses = model_output[:6]
        elif isinstance(model_output, tuple):
            predictions = model_output[0]
            ligs_keypts = recs_keypts = rotations = translations = geom_losses = None
        else:
            predictions = model_output
            ligs_keypts = recs_keypts = rotations = translations = geom_losses = None
        
        # Calculate confidence scores
        confidence_scores = []
        if (ligs_keypts is not None and recs_keypts is not None and 
            rotations is not None and translations is not None):
            
            batch_size = len(predictions)
            for i in range(batch_size):
                try:
                    lig_kpt = ligs_keypts[i]
                    rec_kpt = recs_keypts[i]
                    rotation = rotations[i]
                    translation = translations[i]
                    
                    # Ensure tensors are 2D
                    if len(lig_kpt.shape) == 1:
                        lig_kpt = lig_kpt.unsqueeze(0)
                    if len(rec_kpt.shape) == 1:
                        rec_kpt = rec_kpt.unsqueeze(0)
                    
                    # Calculate keypoint alignment loss
                    transformed_lig = (rotation @ lig_kpt.t()).t() + translation
                    keypoint_alignment_loss = torch.mean(torch.norm(transformed_lig - rec_kpt, dim=1))
                    
                    # Handle geometric loss
                    geom_loss_val = float(geom_losses) if isinstance(geom_losses, (int, float)) else 0.0
                    
                    # Combine losses as confidence indicator (higher = more confident)
                    confidence = -float(keypoint_alignment_loss.cpu()) - geom_loss_val
                    confidence_scores.append(confidence)
                        
                except Exception as e:
                    print(f"Warning: Could not calculate confidence for ligand {i}: {e}")
                    confidence_scores.append(0.0)
        else:
            # Fallback: use dummy confidence scores
            confidence_scores = [0.0] * len(predictions)
        
        out_ligs = ligs
        out_lig_coords = lig_coords
        names = []
        for i, lig in enumerate(ligs):
            # Try to get the ligand name, fallback to index-based naming
            if lig.HasProp("_Name") and lig.GetProp("_Name").strip():
                name = lig.GetProp("_Name")
            else:
                # Use the true_indices for consistent naming
                name = f"{true_indices[i]}"
            names.append(name)
        
        successes = list(zip(true_indices, names, confidence_scores))
        failures = []
        
    except AssertionError:
        # Handle individual processing when batch fails
        lig_graphs, rec_graphs, geometry_graphs = (dgl.unbatch(lig_graphs),
        dgl.unbatch(rec_graphs), dgl.unbatch(geometry_graphs))
        predictions = []
        out_ligs = []
        out_lig_coords = []
        successes = []
        failures = []
        confidence_scores = []
        
        for lig, lig_coord, lig_graph, rec_graph, geometry_graph, true_index in zip(
            ligs, lig_coords, lig_graphs, rec_graphs, geometry_graphs, true_indices):
            try:
                model_output = model(lig_graph, rec_graph, geometry_graph)
                
                if isinstance(model_output, tuple) and len(model_output) >= 6:
                    pred, lig_kpt, rec_kpt, rotation, translation, geom_loss = model_output[:6]
                    prediction = pred[0]
                    
                    try:
                        # Calculate confidence for individual ligand
                        if len(lig_kpt[0].shape) == 1:
                            lig_kpt_tensor = lig_kpt[0].unsqueeze(0)
                        else:
                            lig_kpt_tensor = lig_kpt[0]
                        
                        if len(rec_kpt[0].shape) == 1:
                            rec_kpt_tensor = rec_kpt[0].unsqueeze(0)
                        else:
                            rec_kpt_tensor = rec_kpt[0]
                        
                        transformed_lig = (rotation[0] @ lig_kpt_tensor.t()).t() + translation[0]
                        keypoint_alignment_loss = torch.mean(torch.norm(transformed_lig - rec_kpt_tensor, dim=1))
                        geom_loss_val = float(geom_loss) if isinstance(geom_loss, (int, float, torch.Tensor)) else 0.0
                        confidence = -float(keypoint_alignment_loss.cpu()) - geom_loss_val
                    except Exception as e:
                        print(f"Warning: Could not calculate individual confidence: {e}")
                        confidence = 0.0
                else:
                    prediction = model_output[0][0] if isinstance(model_output, tuple) else model_output[0]
                    confidence = 0.0
                
                out_ligs.append(lig)
                out_lig_coords.append(lig_coord)
                predictions.append(prediction)
                confidence_scores.append(confidence)
                
                # Try to get the ligand name, fallback to index-based naming
                if lig.HasProp("_Name") and lig.GetProp("_Name").strip():
                    name = lig.GetProp("_Name")
                else:
                    name = f"{true_index}"
                    
                successes.append((true_index, name, confidence))
                
            except Exception as e:
                # Try to get the ligand name, fallback to index-based naming
                if lig.HasProp("_Name") and lig.GetProp("_Name").strip():
                    name = lig.GetProp("_Name")
                else:
                    name = f"{true_index}"
                    
                failures.append((true_index, name))
                print(f"Error processing {name}: {e}")
                
    assert len(predictions) == len(out_ligs) == len(confidence_scores)
    return out_ligs, out_lig_coords, predictions, successes, failures, confidence_scores

def run_corrections(lig, lig_coord, ligs_coords_pred_untuned):
    input_coords = lig_coord.detach().cpu()
    prediction = ligs_coords_pred_untuned.detach().cpu()
    lig_input = deepcopy(lig)
    conf = lig_input.GetConformer()
    for i in range(lig_input.GetNumAtoms()):
        x, y, z = input_coords.numpy()[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    lig_equibind = deepcopy(lig)
    conf = lig_equibind.GetConformer()
    for i in range(lig_equibind.GetNumAtoms()):
        x, y, z = prediction.numpy()[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    coords_pred = lig_equibind.GetConformer().GetPositions()

    Z_pt_cloud = coords_pred
    rotable_bonds = get_torsions([lig_input])
    new_dihedrals = np.zeros(len(rotable_bonds))
    for idx, r in enumerate(rotable_bonds):
        new_dihedrals[idx] = get_dihedral_vonMises(lig_input, lig_input.GetConformer(), r, Z_pt_cloud)
    optimized_mol = apply_changes(lig_input, new_dihedrals, rotable_bonds)
    optimized_conf = optimized_mol.GetConformer()
    coords_pred_optimized = optimized_conf.GetPositions()
    
    try:
        R, t = rigid_transform_Kabsch_3D(coords_pred_optimized.T, coords_pred.T)
        coords_pred_optimized = (R @ (coords_pred_optimized).T).T + t.squeeze()
        for i in range(optimized_mol.GetNumAtoms()):
            x, y, z = coords_pred_optimized[i]
            optimized_conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    except np.linalg.LinAlgError:
        # If SVD fails, skip the rigid transformation and use the optimized coordinates as-is
        print("Warning: SVD failed in rigid transformation, using optimized coordinates without alignment")
        pass
    
    return optimized_mol

# def write_while_inferring(dataloader, model, args):
    
#     full_output_path = os.path.join(args.output_directory, "output.sdf")
#     full_failed_path = os.path.join(args.output_directory, "failed.txt")
#     full_success_path = os.path.join(args.output_directory, "success.txt")

#     w_or_a = "a" if args.skip_in_output else "w"
#     with torch.no_grad(), open(full_output_path, w_or_a) as file, open(
#         full_failed_path, "a") as failed_file, open(full_success_path, w_or_a) as success_file:
#         with Chem.SDWriter(file) as writer:
#             i = 0
#             total_ligs = len(dataloader.dataset)
#             for batch in dataloader:
#                 i += args.batch_size
#                 print(f"Entering batch ending in index {min(i, total_ligs)}/{len(dataloader.dataset)}")
#                 ligs, lig_coords, lig_graphs, rec_graphs, geometry_graphs, true_indices, failed_in_batch = batch
#                 for failure in failed_in_batch:
#                     if failure[1] == "Skipped":
#                         continue
#                     failed_file.write(f"{failure[0]} {failure[1]}")
#                     failed_file.write("\n")
#                 if ligs is None:
#                     continue
#                 lig_graphs = lig_graphs.to(args.device)
#                 rec_graphs = rec_graphs.to(args.device)
#                 geometry_graphs = geometry_graphs.to(args.device)
                
                
#                 out_ligs, out_lig_coords, predictions, successes, failures = run_batch(model, ligs, lig_coords,
#                                                                                        lig_graphs, rec_graphs,
#                                                                                        geometry_graphs, true_indices)
                
#                 # Apply corrections with error handling for SVD convergence issues
#                 opt_mols = []
#                 correction_failures = []
#                 for i, (lig, lig_coord, prediction) in enumerate(zip(out_ligs, out_lig_coords, predictions)):
#                     try:
#                         opt_mol = run_corrections(lig, lig_coord, prediction)
#                         opt_mols.append(opt_mol)
#                     except np.linalg.LinAlgError as e:
#                         print(f"Warning: SVD did not converge for molecule {lig.GetProp('_Name') if lig.HasProp('_Name') else i}: {e}")
#                         # Use the original molecule without corrections
#                         opt_mols.append(lig)
#                         correction_failures.append(f"molecule_{i}_svd_failed")
#                     except Exception as e:
#                         print(f"Warning: Error in corrections for molecule {lig.GetProp('_Name') if lig.HasProp('_Name') else i}: {e}")
#                         # Use the original molecule without corrections
#                         opt_mols.append(lig)
#                         correction_failures.append(f"molecule_{i}_correction_failed")
                
#                 for mol, success in zip(opt_mols, successes):
#                     writer.write(mol)
#                     success_file.write(f"{success[0]} {success[1]}")
#                     success_file.write("\n")
#                     # print(f"written {mol.GetProp('_Name')} to output")
#                 for failure in failures:
#                     failed_file.write(f"{failure[0]} {failure[1]}")
#                     failed_file.write("\n")

# Support for extraction and writing of confidence scores
def write_while_inferring(dataloader, model, args):
    
    full_output_path = os.path.join(args.output_directory, "output.sdf")
    full_failed_path = os.path.join(args.output_directory, "failed.txt")
    full_success_path = os.path.join(args.output_directory, "success.txt")
    confidence_scores_path = os.path.join(args.output_directory, "confidence_scores.csv")

    w_or_a = "a" if args.skip_in_output else "w"
    
    # Initialize confidence scores CSV
    import pandas as pd
    confidence_data = []
    
    with torch.no_grad(), open(full_output_path, w_or_a) as file, open(
        full_failed_path, "a") as failed_file, open(full_success_path, w_or_a) as success_file:
        with Chem.SDWriter(file) as writer:
            i = 0
            total_ligs = len(dataloader.dataset)
            for batch in dataloader:
                i += args.batch_size
                print(f"Entering batch ending in index {min(i, total_ligs)}/{len(dataloader.dataset)}")
                ligs, lig_coords, lig_graphs, rec_graphs, geometry_graphs, true_indices, failed_in_batch = batch
                
                for failure in failed_in_batch:
                    if failure[1] == "Skipped":
                        continue
                    failed_file.write(f"{failure[0]} {failure[1]}")
                    failed_file.write("\n")
                if ligs is None:
                    continue
                    
                lig_graphs = lig_graphs.to(args.device)
                rec_graphs = rec_graphs.to(args.device)
                geometry_graphs = geometry_graphs.to(args.device)
                
                out_ligs, out_lig_coords, predictions, successes, failures, confidence_scores = run_batch(
                    model, ligs, lig_coords, lig_graphs, rec_graphs, geometry_graphs, true_indices)
                
                # Apply corrections with error handling for SVD convergence issues
                opt_mols = []
                correction_failures = []
                for i, (lig, lig_coord, prediction) in enumerate(zip(out_ligs, out_lig_coords, predictions)):
                    try:
                        if args.run_corrections:
                            opt_mol = run_corrections(lig, lig_coord, prediction)
                        else:
                            opt_mol = lig
                        opt_mols.append(opt_mol)
                    except np.linalg.LinAlgError as e:
                        print(f"Warning: SVD did not converge for molecule {lig.GetProp('_Name') if lig.HasProp('_Name') else i}: {e}")
                        opt_mols.append(lig)
                        correction_failures.append(f"molecule_{i}_svd_failed")
                    except Exception as e:
                        print(f"Warning: Error in corrections for molecule {lig.GetProp('_Name') if lig.HasProp('_Name') else i}: {e}")
                        opt_mols.append(lig)
                        correction_failures.append(f"molecule_{i}_correction_failed")
                
                # Write molecules and collect confidence data
                for mol, success, conf_score in zip(opt_mols, successes, confidence_scores):
                    writer.write(mol)
                    success_file.write(f"{success[0]} {success[1]}")
                    success_file.write("\n")
                    
                    # Store confidence data
                    ligand_name = success[1]
                    # Ensure .sdf extension
                    if not ligand_name.lower().endswith('.sdf'):
                        filename = f"{ligand_name}.sdf"
                    else:
                        filename = ligand_name
                    
                    confidence_data.append({
                        'filename': filename,
                        'confidence_score': conf_score,
                        'status': 'success'
                    })
                
                for failure in failures:
                    failed_file.write(f"{failure[0]} {failure[1]}")
                    failed_file.write("\n")
                    
                    # Store failed ligand data
                    ligand_name = failure[1]
                    # Ensure .sdf extension
                    if not ligand_name.lower().endswith('.sdf'):
                        filename = f"{ligand_name}.sdf"
                    else:
                        filename = ligand_name
                    
                    confidence_data.append({
                        'filename': filename,
                        'confidence_score': None,
                        'status': 'failed'
                    })
    
    # Save confidence scores to CSV
    if confidence_data:
        df = pd.DataFrame(confidence_data)
        df.to_csv(confidence_scores_path, index=False)
        print(f"Confidence scores saved to: {confidence_scores_path}")

def main(arglist = None):
    args, cmdline_args = parse_arguments(arglist)
    
    args = get_default_args(args, cmdline_args)
    assert args.output_directory, "An output directory should be specified"
    assert args.ligands_sdf, "No ligand sdf specified"
    assert args.rec_pdb, "No protein specified"
    seed_all(args.seed)
    
    os.makedirs(args.output_directory, exist_ok = True)

    success_path = os.path.join(args.output_directory, "success.txt")
    failed_path = os.path.join(args.output_directory, "failed.txt")
    if os.path.exists(success_path) and os.path.exists(failed_path) and args.skip_in_output:
        with open(success_path) as successes, open(failed_path) as failures:
            previous_work = successes.readlines()
            previous_work += failures.readlines()
        previous_work = set(map(lambda tup: int(tup.split(" ")[0]), previous_work))
        print(f"Found {len(previous_work)} previously calculated ligands")
    else:
        previous_work = None
    
        
    rec_graph, model = load_rec_and_model(args)
    if args.lig_slice is not None:
        lig_slice = tuple(map(int, args.lig_slice.split(",")))
    else:
        lig_slice = None
    
    lig_data = multiple_ligands.Ligands(args.ligands_sdf, rec_graph, args, slice = lig_slice, skips = previous_work, lazy = args.lazy_dataload)
    lig_loader = DataLoader(lig_data, batch_size = args.batch_size, collate_fn = lig_data.collate, num_workers = args.n_workers_data_load)

    full_failed_path = os.path.join(args.output_directory, "failed.txt")
    with open(full_failed_path, "a" if args.skip_in_output else "w") as failed_file:
        for failure in lig_data.failed_ligs:
            failed_file.write(f"{failure[0]} {failure[1]}")
            failed_file.write("\n")
    
    write_while_inferring(lig_loader, model, args)

if __name__ == '__main__':
    main()