# Tissuecyte_Annotation
App for Tissuecyte Annotation. Annotations are done in a 10 Micron affine aligned space and then warped to 25 Micron CCF space.

# Environment Setup
First, create a new environment, and then activate the environment (app is currently using python 3.9.12):

```
conda create -n tissuecyte_annotation python=3.9.12
conda activate tissuecyte_annotation
```

After doing this, clone the repo:

```
git clone https://github.com/arjunsridhar12345/Tissuecyte_Annotation.git
```

And then `cd` to the directory where the repo was cloned and install the package and its dependencies

```
pip install -e .
```
Then `cd Software` to begin using the app

# Using the App
### Getting the 10 Micron Volume
First, the 10 micron volume needs to be extracted before proceeding to annotating. To do this, run the following command giving the mouse id that needs to be annotated (example shown below)

```
python resample_images.py --user <allen_institute_username> --password <allen_institute_password> --mouseID <mouseID>
python resample_images.py --user arjun.sridhar --password arjun's password --mouseID 608671
```

Note: Access to HPC/Linux is a requirement as this script calls an executable on the network. To get access to HPC, go to ServiceNow and submit a ticket.
Running this process can take some time (around 20-30 mins usually). The output will be the resampled red, green, and blue images and will be saved to the following directory on the network

`/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID`

### Running the Annotation App
To use the annotation app, run the following command

```
python annotation_app_pyqtgraph_10.py --mouseID <mouseID>
```

The app may take several mintues to load, and once it has loaded, the following screen should be displayed

![image](https://github.com/arjunsridhar12345/Tissuecyte_Annotation/blob/main/images/annotation_app.png)

### Functions of App
  * Point Lock - Disables annotating, deleting, and undoing points. When point lock is ON, all of these features are disabled
  * Annotating - Select the probe and day using the drop downs, and then disable point lock by clicking on the point lock button to set it to OFF. Then, click anywhre on the screen to create an annotation point. Once finished, enable point lock to prevent accidentally annotating unwanted points
  * Annotated points will be automatically saved to a CSV, which has 4 columns (AP, DV, ML, and probe_name). The CSV is stored in the following location
    `/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID/probe_annotations_mouseID.csv`
  * Undo Last Annotation - undos the last annotation made for the selected probe
  * Delete Points for Selected Probe - deletes all of the points for the selected probe on the current slice (only deletes all points on the current slice, not across all slices). A warning message will appear to confirm
  * Delete All Points for Selected Probe - deletes all of the points for the selected probe across all slices. A warning message will appear to confirm
  * Hide Points - Hides the points 
  * Show Points - Shows the points
  * RGB sliders/toggle - modifies the display for the given channels
  * Coronal, Horizontal, Sagittal buttons - change the view depending on which button is pressed
  * To move to different slices arcoss the brain, drag the main slider or use the left/right arrow keys
  * Zoom in/out - use mouse wheel
  * Can move image around my holding the left mouse button and then dragging anywhere on the screen

### Warp Points to CCF
  * Once all points have been annotated, click the Warp Annotations to CCF button to warp the points. Once this process has finished, a popup will appear saying that a new window will open with the warped points. Click OK and then wait for the CCF volume to be loaded. 
  * Then, scroll through the slices to see the warped points. Some of the warped points may appear on different slices, but they should generally be in the same regions of interest when compared with the 10 micron volume.
  * The warped points will be saved to following location
    `/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID/probe_annotations_mouseID_warped.csv`
  * In addition, each of the probes annotated will have a CSV file in the following location
    `/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID/Probe_<probe_day>_annotations_mouseID_warped.csv`

### Reassignment App
This app allows for probes to be reassigned to other probes along with deleting stray annotations. Run the command `python probe_reassignment.py -- mouseID <mouseID> --user <username> --password <password> --implant <implant>` and then wait for some time before seeing a similar display as below. The app shows the probe trajectories in both a 2D and 3D view

![image](https://github.com/arjunsridhar12345/Tissuecyte_Annotation/blob/main/images/reassignment_app_2.png)

The app has the following main components
  * Show/Hide Probes: To show/hide probes, select the checkboxes for the probes that should be hidden, and then click the Update Probe Trajectory Display. Follow the same steps to show the probes by unchecking the desired probes
  * Switch Probes: To switch two probes, select the current probe and day from the corresponding drop down, and then the new probe and day from the corresponding drop down, and then hit the Switch Probes button 
  * Reassign Probes: If a probe needs to be assigned to another one, but you are not sure if it needs to be switched, use this function. For example, if the current probe C1 needs to be reassigned to probe C3, and there is another probe that is C3, using this will change probe C1 to C3, and the original C3 will become the color light grey and have the label origProbe_C3. This allows for reassigning this probe later on, and it provides some history as to what the probe was orignally
  * After reassignments have been done, click the Save Updated Probes button for them to be saved
  
### Generate Images
After the reassignments have been done, open a command prompt terminal and ssh into the hpc using the following command `ssh username@hpc-login` and enter your password. Then follow the steps below

  * To allocate resources, run the command `srun -N1 -c50 -t10:00:00 --mem=250gb -p braintv --pty bash`
  * If allocating takes too long, try changing to --mem=100gb
  * Once resources have been allocated, run this command `cd /allen/scratch/aibstemp/arjun.sridhar/dist/`
  * Finally, run the command `./preprocess_generation --mouseID <mouseID> --useAllProbes yes` if all the probes need to have images/areas generated
  * Otherwise, run the command `./preprocess_generation --mouseID <mouseID>` and only the probes that were reassigned/modified from the previous step will have images generated. This assumes that images have already been generated. If you need to generate all of them again, run the above command
  * If you need to specify a probe, run the following: `./preprocess_generation --mouseID <mouseID> --probe <probe>`
  * The entire process takes around 5 hours to complete if all the probes need to have images generated
  * The same process can be repeated for any other mouse ids by opening another command prompt and following the instructions above
  
### Refinement App
This app allows for the alignment of the 384 channels to the corresponding regions of interest. Run the command `python volume_alignment.py --mouseID <mouseID>`

For dynamic routing experiments, run the command `python volume_alignment_corr.py --mouseID <mouseID> --dr yes`

The probe and the desired metric that will be displayed along with the unit density are shown in the drop downs. The unit density will always be the plot closest to the red probe track.

The app has the following components:
  * Zoom in/out - use the mouse wheel
  * Volume alignment: Zoom in and select a channel (point) on the unit density plot. Then click on a corresponding region on the probe and a yellow horizontal line should connect from the channel to that point on the probe. The plots will then shift/linearly space depending on the alignment made
  * To reposition anchor, click on yellow line and it will turn white as it has been selected. Then use the up/down arrow keys to adjust
  * Remove Alignment: Click on a yellow line. This will turn it white as it has been selected. Then hit the delete button to remove it. Note this works best when zooming out
  * Remove Last Alignment: Hit the left arrow key
  * Reset Metric Plot will clear the alignments and reset both plots to the original positions
  * Toggle Probe will show/hide the line with the ccf regions
  * Toggle Mask will show/hide the CCF mask
  * To change a metric, choose from the drop down
  * To view a different probe, select from the probe drop down
  * When satisfied with alignments, click the Save channel alignments button and a csv will be saved with the ccf coordinate and region for each channel in the following directory `/allen/programs/mindscope/workgroups/np-behavior/tissuecyte/mouseID/Probe_<probe_day>_channels_mouseID_warped.csv`


