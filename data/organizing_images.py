
import math, random
import os, shutil
import matplotlib.pyplot as plt
from PIL import Image

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
from rdkit import Chem

random_mols = random.randint(1, 200)
print(random_mols)

# Config:
images_dir = './output_AURKB'
write_dir = './OUTPUT_AURKB_pics_only'

if not os.path.exists(write_dir):
    os.makedirs(write_dir, exist_ok=True)

images_list = os.listdir(images_dir)
images_list.sort(key=lambda f: int(f.split('_')[-1]))
for folder in images_list:
    print(folder)
    index = folder.split('_')[-1]
    print(index)
    # if index in [58, 59]: next
    for struct in os.listdir(os.path.join(images_dir, folder)):
        if struct.endswith('.sdf'):
            # if index == 59: break
            img_size = (200, 200)
            supplier = Chem.SDMolSupplier(f'{images_dir}/{folder}/{struct}')
            for mol in supplier:
                try:
                    AllChem.Compute2DCoords(mol)
                except:
                    print('Cannot process coordinates')
                    continue
                # mol_id = mol.GetProp('_Name')
                d = rdMolDraw2D.MolDraw2DCairo(*img_size)
                d.DrawMolecule(mol)
                d.FinishDrawing()
                d.WriteDrawingText(f'{index}.png')
                shutil.move(f'{index}.png', write_dir)


# [os.remove(file) for file in os.listdir(images_dir) if file.endswith('.sdf')]

images_dir = write_dir #'/home/luna/Documents/Coding/EquiBind-XAI/onlypics_1k/'
result_grid_location = '.'
result_figsize_resolution = 100 #40 # 1 = 100px
fontsize = 75
images_list = os.listdir(images_dir)
images_list.sort(key=lambda f: int(f.split('.')[0]))
images_count = len(images_list)
print('Images: ', images_list)
print('Images count: ', images_count)

# Calculate the grid size:
grid_size = math.ceil(math.sqrt(25))

# Create plt plot:
# fig, axes = plt.subplots(grid_size, grid_size, figsize=(result_figsize_resolution, result_figsize_resolution))


def plotting(i, images_count):

    plt.clf()
    fig, axes = plt.subplots(grid_size, grid_size, figsize=(result_figsize_resolution, result_figsize_resolution))
    current_file_number = 0
    current_image_number = i
    for image_filename in images_list[i:i+25]:
        x_position = current_file_number % grid_size
        y_position = current_file_number // grid_size

        # plt_image = plt.imread(images_dir + '/' + images_list[current_file_number])
        im_open = Image.open(images_dir + '/' + images_list[current_image_number])
        # width, height = im_open.size 
        # left = 0
        # top = 0
        # right = width / 3
        # bottom = height

        # # takes in (x, y, width, height)
        # im1 = im_open.crop((left, top, right, bottom))
        im1 = im_open
        axes[x_position, y_position].imshow(im1)
        axes[x_position, y_position].set_title(image_filename, fontdict={'fontsize':fontsize})
        axes[x_position, y_position].axis('off')
        print((current_image_number + 1), '/', images_count, ': ', image_filename)

        current_file_number += 1
        current_image_number += 1

    plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.9)
    plt.savefig(os.path.join(result_grid_location, f'grid{i}_{i+25}.png'))

    return


tmp_images_count = 100
for i in range(0, tmp_images_count, 25):
    plotting(i, images_count)