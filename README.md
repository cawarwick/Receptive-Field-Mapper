# Receptive-Field-Mapper
Script which takes the image files from ThorImage and extracts the stimulus events (from ThorSync) to generate a Receptive field map to assist in mapping RFs

After imaging has concluded, move the corresponding Thorsync file into the same folder as the TIFF and Xml files. This is just to simplify the organization as the thorsync files are not named in a derivative fashion of the imaging data. 

<img width="593" height="313" alt="image" src="https://github.com/user-attachments/assets/d32882fe-2bb2-44d5-ab02-618d2b20a156" />

Once the 3 relevant files are all together open the script and make sure it's looking in the right folder. You can copy it directly without replacing the backslashes

<img width="424" height="101" alt="image" src="https://github.com/user-attachments/assets/b9d70a99-4c22-49ea-b561-736402cd19f9" />

Specify the pattern of stimulation which corresponds with the XY program. Not currently able to read PI's Mikromove so just make sure you're using the same pattern. 
Each of these patterns assumes your 0,0 starting point is at top left of the skin (0,0) and are ONLY moving right and down (+X,-Y) 

<img width="638" height="262" alt="image" src="https://github.com/user-attachments/assets/c6d6747b-d9d7-475d-b053-00f9b94af2fd" />

Note that this is expecting 2 channel, multi-Z, time course data with a mechanical stimulator so it doesn't translate to other stimulus channels in thorsync.

How it broadly works
1. extract_TS_metadata() ##extract metadata of image files from the experiment.xml file in order to correctly re-assemble the stacks
2. Load_OMEtiff(Meta_Data) #use metadata to load the OME Tiffs and rearrange the images and save a MaxZ for use in RF analysis
3. Load_ThorSync() ##Load thorSync file for finding stimulus peaks and save the mechanical stimulator's Force channel.
4. Find_Peaks() ##find the start and stop frames of the force based on peak derivative 
5. analyze_stack() #Make image montage of each stimulus (average of frames during stimulus+10 frames) and arrange into RF map.
6. analyze_masked_trace() ##make rough ROIs and extract traces to make heat map of receptive fields. 

Example of RF_Montage for a 3x3 grid area along with max intensity projection to verify the FOV in the SDH is aligned with RF correctly.
<img width="3240" height="2632" alt="RF_montage" src="https://github.com/user-attachments/assets/bfab33e0-7ef1-4277-aa24-50f1a1e9fe20" />

ROI selection. Blob based selection with gaussian smoothing.
<img width="1800" height="855" alt="RF_mask" src="https://github.com/user-attachments/assets/5703698b-4809-42f8-9c21-ae3a2f4f4756" />


Example of RF_ROI_Heatmap. Quantification of all blobs/ROIs found. 
<img width="1500" height="1350" alt="RF_roi_heatmap" src="https://github.com/user-attachments/assets/d48440e1-5aef-4839-bb7c-ddf2983b0c8d" />
