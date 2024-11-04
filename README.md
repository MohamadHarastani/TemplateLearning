# TemplateLearning
Template Learning: Deep Learning with Domain Randomization for Particle Picking in Cryo-Electron Tomography


## Software
### pre-requisites 
- IMOD (we tested with v 4.11.15, but we believe that any version should be compatible)
- Conda (miniconda should be fine)
- Cuda (we tested with 10.2 and 11.4, if you manage to install Eman2 and Parakeet you should be fine)
### Installation
Make sure you have IMOD and Cuda executables in your path
```
which newstack tilt trimvol
which nvcc
```
clone and install the package
```
git clone https://github.com/MohamadHarastani/TemplateLearning.git
conda env create -f environment.yml
```
The following command usually installs everything. If you face any problem, install step by step as follows
```
conda update --all
conda install mamba -c conda-forge
mamba create -y -n TemplateLearning python=3.9 eman-dev==2.99.47 -c cryoem -c conda-forge
conda activate TemplateLearning
pip install ProDy==2.4.1 mrcfile==1.4.3 scikit-image==0.21.0 pyfftw==0.13.1 python-parakeet==0.4.5 h5py==3.8.0
```
## Short tutorial (a detailed tutorial will be provided soon)
- Replace the templates (in PDB or CIF format) in the directory "templates" with templates for your molecule. The current templates are nucleosomes.
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
- Create dense distribution of your volume. This algorithm "Tetris" generates dense distributions placing your templates and distractors very close to each other
```
python step4_CreateTetris.py  # with default values should take ~ 1 hour 15 minutes to create 48 simulated volumes (on a basic CPU)
```
Validation: see the output of step 4 by opening some tetris volumes in the folder "tetrises", make sure you can see your template in between the distractors. Also, make sure the tetris looks dense!
- Start the simulation of the data using the physics simulator (Parakeet)
```
python step5_SimulateData.py >> log.txt  # Depending on the speed of your GPU, on a single GPU, this step will take around 17 hours with default values. If you have multiple GPUs, the time will be reduce linearly
```
Validation: during the simulation, you will get an estimate on how much time remaining after the first iteration. After one iteration is done, you can open the simulated tiltseries with IMOD. The image called 'optics.h5' is very handy to visualize, as it is not too noisy. If it looks empty, you may need to reload the image to see it (IMOD>Edit>Image>Reload>Calc&Apply). It should be simular to the tetris (but tilt series version with defocus) but with a sampling of 1 A/pix. If it doesn't correspond to a tetris, then probably there is a mistake in setting the sizes. Otherwise, sit back and relax!
- Extract what you need to train your model
```
python BinReorderReconstruct.py  # should take a few minutes
```
Validation: look at the folder called "results". It should have tomograms, segmentation maps and coordinates that can be used to train deeplearning models.
- To use deepfinder in the easiest way, use the Scipion workflow (DeepFinder_Scipion.json) attached to the project. Open Scipion, Project -> Import workflow and just populate the protocols with your simulated and expiremental data!
- If you wish to use DeepFinder outside of Scipion, install it from [cryoet-deepfinder]([URL](https://github.com/deep-finder/cryoet-deepfinder)) then run the following step:
```
python Step7_prepareDataForDeepFinder.py  # should be instant
conda activate dfinder # assuming this is the name for the environement
train -p DeepFinder/params_train.xml 
```

Enjoy, and for any question open a ticket!
