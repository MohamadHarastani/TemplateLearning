import xml.etree.ElementTree as ET
import os
import random

# Function to convert text files to XML
def convert_txt_files_to_xml(input_files, output_file):
    objlist = ET.Element("objlist")
    for input_file in input_files:
        tomoID = int(os.path.basename(input_file).split('.')[0])
        with open(input_file, 'r') as f:
            for line in f:
                x, y, z = line.strip().split()
                obj = ET.SubElement(objlist, "object")
                obj.set("tomo_idx", str(tomoID))
                obj.set("class_label", "1")
                obj.set("x", x)
                obj.set("y", y)
                obj.set("z", z)
    tree = ET.ElementTree(objlist)
    tree.write(output_file)

# Paths
output_dir = "DeepFinder/"
os.makedirs(output_dir, exist_ok=True)

NumberOfTomograms = 48

# Random validation ID selection and file list creation
validation_ID = random.randint(0, NumberOfTomograms - 1)
training_files = [f"results/coordinates/{i}.txt" for i in range(NumberOfTomograms) if i != validation_ID]
convert_txt_files_to_xml(training_files, os.path.join(output_dir, "objl_train.xml"))

validation_file = [f"results/coordinates/{validation_ID}.txt"]
convert_txt_files_to_xml(validation_file, os.path.join(output_dir, "objl_valid.xml"))

# Create the params_train.xml
params_train = ET.Element("paramsTrain")
ET.SubElement(params_train, "path_out", path=output_dir)

# Add tomograms paths
path_tomo = ET.SubElement(params_train, "path_tomo")
for i in range(NumberOfTomograms):
    ET.SubElement(path_tomo, f"tomo{i}", path=f"results/tomograms/{i}.mrc")

# Add segmentations paths
path_target = ET.SubElement(params_train, "path_target")
for i in range(NumberOfTomograms):
    ET.SubElement(path_target, f"target{i}", path=f"results/segmentations/{i}.mrc")

# Add training and validation object list paths
ET.SubElement(params_train, "path_objl_train", path=os.path.join(output_dir, "objl_train.xml"))
ET.SubElement(params_train, "path_objl_valid", path=os.path.join(output_dir, "objl_valid.xml"))

# Additional training parameters
ET.SubElement(params_train, "number_of_classes", n="2")
ET.SubElement(params_train, "patch_size", n="48")
ET.SubElement(params_train, "batch_size", n="25")
ET.SubElement(params_train, "number_of_epochs", n="100")
ET.SubElement(params_train, "steps_per_epoch", n="100")
ET.SubElement(params_train, "steps_per_validation", n="10")
ET.SubElement(params_train, "flag_direct_read", flag="False")
ET.SubElement(params_train, "flag_bootstrap", flag="True")
ET.SubElement(params_train, "random_shift", shift="13")

# Write the params_train.xml file
tree = ET.ElementTree(params_train)
tree.write("DeepFinder/params_train.xml")

print("params_train.xml and objlist files created successfully.")


