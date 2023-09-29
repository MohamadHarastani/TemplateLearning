# TemplateLearning
Template Learning: Cryo-ET Particle-picking by Training Deep Learning on Templates via Domain Randomization


## Software
### pre-requisites 
- Imod (we tested with v 4.11.15, but we believe that any version should be compatible)
- Conda (miniconda should be fine)
- Cuda (we tested with 10.2 and 11.4, if you manage to install Eman2 and Parakeet you should be fine)
### Installation
Make sure you have imod installed and in the path
```
which newstack tilt trimvol
```
Make sure you have Cuda in your path
```
which nvcc
```
clone and install the package
```
git clone https://github.com/MohamadHarastani/TemplateLearning.git
conda env create -f environment.yml
```
If the above command doesn't work, try installing the following and fix the conflicts
```
conda update --all
conda install mamba -c conda-forge
mamba create -y -n TemplateLearning python=3.9 eman-dev==2.99.47 -c cryoem -c conda-forge
conda activate TemplateLearning
pip install ProDy==2.4.1 mrcfile==1.4.3 scikit-image==0.21.0 pyfftw==0.13.1 python-parakeet==0.4.5 h5py==3.8.0
```
## Short tutorial (see the detailed tutorial in case of any doubt)
- Replace the templates (in PDB or CIF format) in the directory "templates" with templates for your molecule
- If you wish to simulate additional conformational variability (recommended), use steps 1 and 2. Otherwise, skip to step 3
```
conda activate TemplateLearning
python step1_NmaCalculate.py  # with default values, should take a few seconds for small templates (~200 kDa), up to a few hours for huge templates (~2 MDa)
python step2_DeformTemplates.py # with default values, should not take more than a few minutes
```
Validation: see the output of step 2 by opening some simulated structures in the folder "templates"
- Create volumes from your templates and distractors using step 3
```
python step3_CreateVolumes.py  # with default values should take ~ 10 minutes
```
Validation: see the output of step 3 by opening some volumes in the folder "volumes", especially those for your templates.
- Create dense distribution of your volume. This algorithm "snowball" generates dense distributions placing your templates and distractors very close to each other
```
python step4_CreateSnowballs.py  # with default values should take ~ 1 hour 15 minutes (on a basic CPU)
```
Validation: see the output of step 4 by opening some snowball volumes in the folder "snowballs", make sure you can see your template in between the distractors. Also, make sure the snowball looks dense!
- Start the simulation of the data using the physics simulator
```
python step5_SimulateData.py >> log.txt  # sit back and relax! Depending on the speed of your GPU, this step will take around 17 hours with default values.
```
Validation: during the simulation, you will get an estimate on how much time remaining after the first iteration. After one iteration is done (by default parakeet/32), you can open the simulated tiltseries with imod. The image called 'optics.h5' is very handy to visualize, as it is not too noisy. If it looks empty, you may need to reload the image to see it (Imod>Edit>Image>Reload>Calc&Apply). It should be simular to the snowball but with a sampling of 1 A/pix. If it doesn't correspond to a snowball, them probably you have a mistake in setting the sizes. Otherwise, keep the simulation going
- Extract what you need to train your model
```
python BinReorderReconstruct.py  # should take a few minutes
```
Validation: look at the folder called "results". It should have tomograms, segmentation maps and coordinates that can be used to train deeplearning models.
- To use deepfinder in the easiest way, use the Scipion workflow (DeepFinder_Scipion.json) attached to the project. Open Scipion, Project -> Import workflow
