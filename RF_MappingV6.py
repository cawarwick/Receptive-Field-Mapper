import xml.etree.ElementTree as ET
from pathlib import Path
import tifffile as tif
import numpy as np
import os
import math
import h5py
import csv
import os
from scipy.signal import resample
from scipy import signal
import scipy.signal as signal
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy import ndimage as ndi
import subprocess

# ---------------- USER PARAMETERS ---------------- #

Imaging_Folder = r"E:\PrRFMap\#550_009"

# Physical (X, Y) location of each stimulus, in presentation order (1st event -> 1st entry, etc.)
Pattern_5mmx5mm_1mm_spacing_6x6 = [
    (0, 0), (3, -3), (0, -1), (3, -4),
    (4, 0), (2, -3), (5, 0), (1, -1),
    (2, -5), (3, -1), (0, -4), (4, -5),
    (1, -2), (5, -3), (2, 0), (1, -4),
    (4, -1), (0, -2), (3, -5), (5, -1),
    (2, -4), (3, 0), (4, -4), (0, -3),
    (4, -2), (1, -5), (2, -1), (5, -4),
    (1, -3), (5, -2), (1, 0), (4, -3),
    (0, -5), (3, -2), (5, -5), (2, -2),
]
Pattern_5mmx5mm_2_5mm_spacing_3x3 = [
    (0.0,  0.0),    # 1
    (5.0, -2.5),    # 2
    (0.0, -5.0),    # 3
    (2.5,  0.0),    # 4
    (2.5, -5.0),    # 5
    (5.0,  0.0),    # 6
    (0.0, -2.5),    # 7
    (5.0, -5.0),    # 8
    (2.5, -2.5),    # 9
]
Pattern_7mmx7mm_1mm_spacing_8x8 = [
    (0, 0), (4, -4), (0, -1), (5, -4),
    (0, -7), (4, -3), (1, -7), (4, -2),
    (0, -6), (5, -3), (0, -5), (5, -2),
    (2, -7), (3, -1), (7, -5), (2, -2),
    (6, -6), (1, -3), (7, -2), (2, -5),
    (5, 0), (3, -5), (6, 0), (2, -4),
    (7, -1), (1, 0), (4, -5), (2, 0),
    (5, -5), (1, -1), (4, -6), (4, 0),
    (1, -5), (6, -2), (3, -7), (0, -2),
    (6, -3), (1, -6), (4, -1), (5, -7),
    (5, -1), (4, -7), (1, -2), (6, -5),
    (0, -4), (6, -4), (0, -3), (5, -6),
    (2, -1), (7, -4), (1, -4), (7, -3),
    (2, -6), (3, 0), (3, -6), (6, -1),
    (7, -7), (3, -3), (6, -7), (3, -2),
    (7, -6), (2, -3), (7, 0), (3, -4),
]
Pattern_10mmx10mm_2_5mm_spacing_5x5 = [
    (0, 0), (5, -7.5), (2.5, 0), (10, -5),
    (2.5, -10), (7.5, -2.5), (0, -7.5), (5, 0),
    (7.5, -7.5), (0, -5), (7.5, -10), (2.5, -2.5),
    (10, -7.5), (7.5, 0), (2.5, -7.5), (10, -2.5),
    (5, -10), (0, -2.5), (7.5, -5), (0, -10),
    (5, -2.5), (10, -10), (2.5, -5), (10, 0),
    (5, -5),
]
Pattern_3mmx3mm_1mm_spacing_4x4 = [
    (0, 0), (1, -2), (2, 0), (0, -1),
    (1, -3), (3, -2), (3, 0), (2, -2),
    (0, -3), (1, -1), (2, -3), (3, -1),
    (1, 0), (0, -2), (2, -1), (3, -3),
]
Pattern_7_5mmx7_5mm_2_5mm_spacing_4x4 = [
    (0, 0), (5, -2.5), (2.5, -7.5), (0, -2.5),
    (5, 0), (7.5, -5), (7.5, 0), (5, -5),
    (0, -7.5), (2.5, -2.5), (5, -7.5), (0, -5),
    (2.5, 0), (7.5, -2.5), (2.5, -5), (7.5, -7.5),
]

# Copy paste the desired patttern below into Stim_XY to selet the program. 
STIM_XY = Pattern_5mmx5mm_2_5mm_spacing_3x3


# -------------------------------------------------- #

#Main body of code#
 
# -------------------------------------------------- #

##this is just redundant with above, it's just easier to read without the XY coorindatese.
options = {
    "5mmx5mm 1mm spacing 6x6": Pattern_5mmx5mm_1mm_spacing_6x6,
    "5mmx5mm 2.5mm spacing 3x3": Pattern_5mmx5mm_2_5mm_spacing_3x3,
    "7mmx7mm 1mm spacing 8x8": Pattern_7mmx7mm_1mm_spacing_8x8,
    "10mmx10mm 2.5mm spacing 5x5": Pattern_10mmx10mm_2_5mm_spacing_5x5,
    "3mmx3mm 1mm spacing 4x4": Pattern_3mmx3mm_1mm_spacing_4x4,
    "7.5mmx7.5mm 2.5mm spacing 4x4": Pattern_7_5mmx7_5mm_2_5mm_spacing_4x4
}


Imaging_Folder = (Imaging_Folder.replace(os.sep, "/")+"/") 
split= Imaging_Folder.split('/')
num_folders = len(split)
Parent_Folder_Name=split[(num_folders-2)]

def extract_TS_metadata():
    ##this function extracts metadata from the ThorImageLS XML file and returns it as a dictionary for use in re-arrangin the OME TIFFs
    ##The experiment.xml must be in the same folder as the OME TIFFs for this to work.
    location = Imaging_Folder+"Experiment.xml"
    tree = ET.parse(location)
    root = tree.getroot()

    Meta_Data = {
        "steps": None,
        "pixelX": None,
        "pixelY": None,
        "frameRate": None,
        "flybackCycles": None,##This is actually # Z plus # flyback, so total frame cycles
        "Timelapse timepoints": None,
        "pixelSizeUM": None,
        "ChannelEnable Set": None,
    }

    for elem in root.iter():
        tag = elem.tag.split("}")[-1]

        # ----- Standard tag text extraction -----
        if elem.text and elem.text.strip():
            if tag in Meta_Data:
                Meta_Data[tag] = elem.text.strip()

        # ----- Attribute-based extraction -----
        for attr, value in elem.attrib.items():

            # Timelapse timepoints
            if tag == "Timelapse" and attr == "timepoints":
                value = int(value)
                Meta_Data["Timelapse timepoints"] = value
                
            # ChannelEnable Set
            if tag == "ChannelEnable" and attr == "Set":
                Meta_Data["ChannelEnable Set"] = value
                value = int(value)
                if value == 1:##Note I can't figure out how to detect channels based only on the metadata. This needs to be either fixed or removed entirely.
                    Meta_Data["ChannelEnable Set"] = 1
                if value == 3:
                    Meta_Data["ChannelEnable Set"] = 2
                if value == 11:
                    Meta_Data["ChannelEnable Set"] = 3
                else:
                    print("Channels were not detected as 1,2, or 3 double check input data")

            # Generic attribute fallback
            if attr in Meta_Data:
                Meta_Data[attr] = value

    return Meta_Data

def Load_OMEtiff(Meta_Data):
    ###this function loads the OME TIFFs and re-arranges them into the right order for typical analysis 
    # It also saves a green channel maximum intensity projection in a subfolder.
    ## Note that this assumes there are two channels since I couldn't find the parameter in the metadata that tells you how many channels there are. 
    print("\nTIFF Metadata:")
    print("-" * 40)
    for k, v in Meta_Data.items():
        print(f"{k:25}: {v}")
    Z = Meta_Data.get("steps")

    X = int(Meta_Data['pixelX'])   # pixels
    Y = int(Meta_Data['pixelY'])
    Z = int(Meta_Data['steps'])
    T = int(Meta_Data['Timelapse timepoints'])
    C = 2 #Meta_Data['pixelX']

    print("Z steps", Z)
    print("flyback cycles", Meta_Data['flybackCycles'])

    MAX_FILE_SIZE = 3.8 * 1024**3  # 4 GB
    
    # --------------------------------------------------
    # LOAD TIFF (CHECK FOR MULTIPLE REGIONS)
    # --------------------------------------------------
    input_tif = Imaging_Folder+"Image_scan_1_region_0_0.tif" ##this is the current namaing for ThorImageLS

    # Check for second file
    base, ext = os.path.splitext(input_tif)

    second_tif = Imaging_Folder+"Image_scan_1_region_0_1.tif" ##this is the current namaing for ThorImageLS

    if os.path.exists(second_tif):

        print(f"Multiple TIFF files detected:")
        print(f"  {input_tif}")
        print(f"  {second_tif}")

        with tif.TiffFile(input_tif) as tf:
            stack1 = tf.asarray()

        with tif.TiffFile(second_tif) as tf:
            stack2 = tf.asarray()

        # Concatenate along the image-plane dimension
        stack = np.concatenate((stack1, stack2), axis=0)

        print(f"Combined stack shape: {stack.shape}")

    else:

        print("No second TIFF found. Loading single file.")

        with tif.TiffFile(input_tif) as tf:
            stack = tf.asarray()

    #print("TIFF loading complete. Continuing processing...")
    
    # Save the name of the parent folder to save the tiffs with the right name
    split= Imaging_Folder.split('/')
    num_folders = len(split)
    Parent_Folder_Name=split[(num_folders-2)]
    #print(Parent_Folder_Name)
    output_folder = Imaging_Folder + "Rearranged_Output/" 
    os.makedirs(output_folder, exist_ok=True)
    output_name = output_folder + Parent_Folder_Name

    total_planes = stack.shape[0]
    expected_planes = Z * T * C

    if total_planes != expected_planes:
        raise ValueError(
            f"Plane count mismatch: found {total_planes}, expected {expected_planes}"
        )

    #print("Reshaping to XYZTC...")
    # Original assumed order: T, C, Z, Y, X
    stack = stack.reshape(C, T, Z, Y, X)

    # --------------------------------------------------
    # CREATE GREEN CHANNEL MAXIMUM INTENSITY PROJECTION
    # --------------------------------------------------

    #print("Generating green channel maximum projection...")

    #print(stack.shape)

    # Stack dimensions: C, T, Z, Y, X
    green_stack = stack[0, :, :, :, :]

    # Maximum intensity projection across Z
    green_mip = np.max(green_stack, axis=1)

    #print(green_mip.shape)
    
    mip_output = f"{output_name}_green_MIP.tif"

    green_mip_mean = np.mean(green_mip, axis=(1, 2))

    #print(f"Saving {green_mip.shape[0]} individual MIP pages...")

    tif.imwrite(
        mip_output,
        green_mip,                      # shape (T, Y, X)
        imagej=True,
        metadata={"axes": "TYX"},
        photometric="minisblack",
        compression=None,
    )
    try:
        imagej_exe = r"C:\Users\warwickc\Desktop\Fiji new\Fiji.app\ImageJ-win64.exe"  # Update this path to your ImageJ executable
        open_in_imagej(mip_output, imagej_exe) ##optional but make sure to specify the path to ImageJ on the local machine
    except Exception as e:
        print(f"Could not open ImageJ: {e}")
   
    #print(f"Green MIP saved: {mip_output}")

    """
    # --------------------------------------------------
    ##For saving the full XYZTC stacks everything below this is only useful for saving the full re-arranged data. 
    # --------------------------------------------------

    print("Reordering to XYCZT (page order C fastest)...")
    # Desired iteration order:
    # T -> Z -> C
    reordered = []
    for t in range(T):
        for z in range(Z):
            for c in range(C):
                reordered.append(stack[c, t, z])

    reordered = np.asarray(reordered)


    bytes_per_plane = reordered[0].nbytes
    planes_per_file_max = int(MAX_FILE_SIZE // bytes_per_plane)

    # Must be divisible by Z*C
    planes_per_block = int(Z * C) 
    planes_per_file = int((planes_per_file_max // planes_per_block) * planes_per_block)

    if planes_per_file == 0:
        raise ValueError("Single Z×C block exceeds 4 GB")

    n_files = math.ceil(len(reordered) / planes_per_file)

    print(f"Writing {n_files} output files...")

    for i in range(n_files):
        start = int(i * planes_per_file)
        end = int(min((i + 1) * planes_per_file, len(reordered)))

        out_name = f"{output_name}_{i+1:03d}.tif"

        print(f"Writing {out_name} ({end - start} planes)")

        tif.imwrite(
            out_name,
            reordered[start:end],
            compression=None,
            bigtiff=False
        )
    """

    return(green_mip_mean,green_mip)

def open_in_imagej(tiff_path, imagej_exe):
    tiff_path = str(Path(tiff_path).resolve())
    subprocess.Popen([imagej_exe, tiff_path])   # Popen = don't block Python

def downsample(data, original_fs, target_fs):
    ##this function downsamples the thorsync data (typically 30khz) to the target frequency. It uses linear interpolation to resample the data.
    # Generate the original and target time vectors
    original_time = np.arange(0, len(data) / original_fs, 1 / original_fs)
    target_time = np.arange(0, len(data) / original_fs, 1 / target_fs)
    
    # Use numpy's interp function to downsample the data
    downsampled_data = np.interp(target_time, original_time, data.flatten())
    
    return downsampled_data

def Load_ThorSync(target_fs=1000.0, flybackCycles=6):
    ##This opens up the ThorSync HDF5 file and extracts ONLY force and frame counters. Frame counter is the number of frames (Z+flyback)*T rather than number of time points*Z steps. 
    ##This assumes the original frame rate is 30khz and it downsamples to 1khz. It also smooths the force trace to remove small amounts of noise 
    ##If there is a different stimuli, make sure the names match excatly as it matches the names to extract the data. 
    file_path = Imaging_Folder + "Episode_0000.h5"
    
    selected_datasets = [
        #"AI/E_Stim",
        "AI/Force",
        #"AI/M_Length",
        #"AI/Pedal",
        #"AI/Pelt_Temp",
        #"AI/PiezoMonitor",
        #"AI/PockelsMonitor",
        #"DI/2pFrames",
        #"DI/CaptureActive",
        #"DI/PandaFrames",
        "CI/FrameCounter"
    ]

    try:
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_file_name = f"{base_name}_selected_datasets.csv"

        with h5py.File(file_path, 'r') as h5_file:
            # Downsample all selected datasets and organize into a NumPy array
            downsampled_data = []

            for dataset_name in selected_datasets:
                item = h5_file[dataset_name]
                data = item[:]
                original_fs = 30000.0  # Replace with the actual sampling frequency of your data
                downsampled_data.append(downsample(data, original_fs, target_fs))
                
            # Transpose the array to have time in the first column and each dataset in subsequent columns
            downsampled_data = np.column_stack(downsampled_data)

            #convert frames to volumes and rearrange so 'time' is in the first axis
            Frames=downsampled_data[:,1]
            Frames=np.divide(Frames,int(flybackCycles))##this is the number of planes plus the number of flyback frames
            downsampled_data[:,1]=downsampled_data[:,0]
            downsampled_data[:,0]=Frames

            # Smooth the force traces
            data_to_smooth=downsampled_data[:,1]
            smoothed_data = signal.savgol_filter(data_to_smooth,41, 0)##Smoothing is strictly to deal with small amounts of noise which dramatically affect the derivative
            downsampled_data[:,1] = smoothed_data

            # Write the NumPy array to a CSV file with the original file name
            header_row_min = ["Frames","Force"]

            output_file_name=(Imaging_Folder+base_name+".csv")
            np.savetxt(output_file_name, downsampled_data, delimiter=",", header=",".join(header_row_min), comments="")
            #np.savetxt(output_file_name, downsampled_data, delimiter=",", header=",".join(header_row), comments="")
            #print(f"CSV file '{output_file_name}' created.")          
    except Exception as e:
        print(f"Error: {e}")
    return  downsampled_data

def Find_Peaks(Frame_Force, Imaging_Folder, threshold=.2, window_length=11, polyorder=0, scaling_factor=30):
    ##take in the downsampled force trace and find the peaks in the derivative of the force trace. It returns a list of median values, start and stop times for each peak.
    ##Frame_Force contains the time course (in frames T frames rather than flyback) and the continuous force trace.  
    ##threshold is in arbitrary numbers, so if things change calibrations this will need to be differnet
    ## There is a second round of smoothing here since Load_ThorSync already smooths the force trace, but it's working so just a note that this is a second round of smoothing on the force traces.
    data = Frame_Force
    time_series = data[:, 0]
    data_series = data[:, 1]

    # Parameters  for lowpass filtering
    fs = 1000  # Sampling frequency in Hz
    cutoff = 40  # Hz (cutoff frequency for filtering). 20hz started to affect the square wave, so 40-100hz is a good range.
    order = 4  # Filter order

    # Design low-pass Butterworth filter
    nyquist = 0.5 * fs  # Nyquist frequency
    normal_cutoff = cutoff / nyquist
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    filtered_data = signal.filtfilt(b, a, data_series)
    
    # Save the lowpass filtered data to a CSV
    #LP_Filtered_Force = pd.DataFrame({'Time': time_series[:], 'LP_Filt': filtered_data})
    #LP_Force_filename = os.path.join(output_dir, f'LP_Filt{os.path.basename(csv_file)}')
    #P_Filtered_Force.to_csv(LP_Force_filename, index=False)

      # Find peaks in the derivative

    threshold = 0.1
    min_duration = 100 ##this is in frames of input data, so if 1khz sampling rate this is 100ms minimum 
    scaling_factor = 40  # Scaling factor for the median value. e.g. it's to convert the volts to mN of force. ~40 is ok, but not quite linear.

    above_threshold = filtered_data > threshold

    # Find where the data transitions above the threshold
    starts = np.where(np.diff(above_threshold.astype(int)) == 1)[0] + 1
    ends = np.where(np.diff(above_threshold.astype(int)) == -1)[0] + 1

    if above_threshold[0]:
        starts = np.insert(starts, 0, 0)  # Case where we start above the threshold
    if above_threshold[-1]:
        ends = np.append(ends, len(filtered_data))  # Case where we end above the threshold

    results = []
    for start, end in zip(starts, ends):
        if (end - start) >= min_duration:
            median_value = np.median(filtered_data[start:end]) * scaling_factor  # Apply scaling factor
            results.append([median_value, time_series[start], time_series[end]])

    results_array = np.array(results)
    #print("Results array (median, start, end):")
    #print(results_array)
    # Save results to CSV

    # Specify the directory to save the file
    save_directory = Imaging_Folder  # Replace with your desired directory path

    # Full path for the CSV file
    new_name="peaks_"+ Parent_Folder_Name + ".csv"
    csv_file_path = os.path.join(save_directory, new_name)
    
    # Save results to CSV
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Median mN', 'Start', 'End'])  # Write header row
        writer.writerows(results_array)  # Write the data rows

    #print(f"Results saved to '{csv_file_path}'.")
    return results_array

def analyze_events(trace, events, stim_xy=STIM_XY, baseline_frames=10,
                   post_frames=0, stop_inclusive=True, cmap="viridis",
                   plot=True, ax=None, save_prefix=None):
    """

    ----------
    trace : 1D array-like, e.g. shape (8000,) or (1, 8000)
    events : (N, 3) array-like; columns = (intensity, start_frame, stop_frame).
        Events are assumed to be in presentation order, so event i maps to
        stim_xy[i % len(stim_xy)]. If there are repeats of the 9-stimulus set,
        repeats at the same location are averaged on the heat map.
    stim_xy : list of (x, y) physical locations in presentation order.
    baseline_frames : frames before start used for the median baseline.
    post_frames : extra frames after stop included in the max window (default 0).
    stop_inclusive : include the stop frame in the max window.
    save_prefix : str or None. If given (e.g. "output/run1"), writes:
        <prefix>_heatmap.png  annotated heat map figure (300 dpi)
        <prefix>_raster.png   bare raster, 1 pixel per stimulus location
        <prefix>_values.csv   per-event raw max, baseline median, ratio
        <prefix>_grid.csv     the 2D ratio grid
 
    Returns dict: intensity, start, stop, peak, baseline, ratio, xy, grid,
    x_vals, y_vals, fig, ax. grid rows = Y (top = max Y), cols = X ascending.
    """
    trace = np.asarray(trace, dtype=float).ravel()
    events = np.atleast_2d(np.asarray(events, dtype=float))
    if events.shape[1] != 3:
        raise ValueError("events must have 3 columns: (intensity, start, stop)")
 
    n_ev, n_stim = len(events), len(stim_xy)
    intensity = events[:, 0]
    starts = events[:, 1].astype(int)
    stops = events[:, 2].astype(int)
 
    peak = np.full(n_ev, np.nan)
    baseline = np.full(n_ev, np.nan)
 
    for i, (s, e) in enumerate(zip(starts, stops)):
        s = max(s, 0)
        e_end = min(e + (1 if stop_inclusive else 0) + post_frames, len(trace))
        if e_end <= s:
            continue
        peak[i] = np.max(trace[s:e_end])
        b0 = max(s - baseline_frames, 0)
        if s - b0 < baseline_frames:
            print(f"Warning: event {i + 1} has only {s - b0} baseline frames available.")
        if s > b0:
            baseline[i] = np.median(trace[b0:s])
 
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = peak / baseline
 
    xy = np.array([stim_xy[i % n_stim] for i in range(n_ev)], dtype=float)
 
    x_vals = np.unique([p[0] for p in stim_xy])
    y_vals = np.unique([p[1] for p in stim_xy])[::-1]
    sums = np.zeros((len(y_vals), len(x_vals)))
    counts = np.zeros_like(sums)
    for (x, y), r in zip(xy, ratio):
        if np.isnan(r):
            continue
        xi = np.where(x_vals == x)[0][0]
        yi = np.where(y_vals == y)[0][0]
        sums[yi, xi] += r
        counts[yi, xi] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        grid = np.where(counts > 0, sums / counts, np.nan)
 
    fig = None
    if plot or save_prefix:
        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 4.5))
        else:
            fig = ax.figure
        im = ax.imshow(grid, cmap=cmap, aspect="equal")
        ax.set_xticks(range(len(x_vals)))
        ax.set_xticklabels([f"{v + 0.0:g}" for v in x_vals])
        ax.set_yticks(range(len(y_vals)))
        ax.set_yticklabels([f"{v + 0.0:g}" for v in y_vals])
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        fig.colorbar(im, ax=ax, label="Peak / baseline median")
        ax.set_title("Response map")
        fig.tight_layout()
 
    if save_prefix:
        d = os.path.dirname(save_prefix)
        if d:
            os.makedirs(d, exist_ok=True)
        fig.savefig(f"{save_prefix}_heatmap.png", dpi=300)
        #plt.imsave(f"{save_prefix}_raster.png", grid, cmap=cmap)
        table = np.column_stack([np.arange(1, n_ev + 1), intensity, starts, stops,
                                 xy[:, 0], xy[:, 1], peak, baseline, ratio])
        np.savetxt(f"{save_prefix}_values.csv", table, delimiter=",", comments="",
                   fmt=["%d", "%g", "%d", "%d", "%g", "%g", "%.6g", "%.6g", "%.6g"],
                   header="event,intensity,start,stop,x,y,max_value,baseline_median,ratio")
        #np.savetxt(f"{save_prefix}_grid.csv", grid, delimiter=",", fmt="%.6g")
        if not plot:
            plt.close(fig)
            fig = None
 
    return {"intensity": intensity, "start": starts, "stop": stops,
            "peak": peak, "baseline": baseline, "ratio": ratio,
            "xy": xy, "grid": grid, "x_vals": x_vals, "y_vals": y_vals,
            "fig": fig, "ax": ax}
 
def analyze_stack(stack, events, stim_xy=STIM_XY, baseline_frames=30,
                  post_frames=10, stop_inclusive=True, stat="mean", 
                  time_axis=0, cmap="magma", blur_sigma=1,
                  symmetric_clim=False, saturation=0.5, vmin=None, vmax=None,
                  mip_saturation=1.0, mip_vmin=0.0, mip_vmax=None,
                  plot=True, save_prefix=None, Savetiff=False):
    """
    #This function takes the discrete start and stop times of the mechanical stimulation (force) and extracts 
    #1. the stimulus response of each stimulus normalized to baseline and maps it into a RF map
    #2. the outputs are a heatmap of each stimulus as well as a montage so you can visualize where the cells are in the FOV
    #3. Also worth noting that the function attempts to find ROIs to improve SNR of the events so if your cells are blurry/dense tweak the ROI method
    Parameters
    
    Image version of analyze_events. For each event:
        delta_image = stimulus image - baseline image
    where
        stimulus image = per-pixel `stat` ("mean" or "max") over frames
                         start .. stop + post_frames
        baseline image = per-pixel mean of the `baseline_frames` frames before start
    A 2D Gaussian blur (blur_sigma) is applied to the images before subtraction.
    ("max" blurs every frame in the window first, then takes the max over time,
    so single-pixel noise spikes are not picked up by the max.)
 
    The per-event images are tiled into a montage using the same stimulus ->
    (x, y) layout as the heat map (X ascending left to right, max Y on the top
    row). Repeats at the same location are averaged. A maximum intensity
    projection (MIP) across all locations is drawn below the montage.
 
    Parameters
    ----------
    stack : ndarray, shape (T, H, W) (or set time_axis if time is elsewhere)
    events : (N, 3) array; columns = (intensity, start_frame, stop_frame).
    baseline_frames : frames before start for the baseline image (default 30).
    post_frames : extra frames after stop included in the stimulus window
        (default 10; try 5 to 10).
    stat : "mean" (default) or "max" over the stimulus window.
    blur_sigma : Gaussian sigma in pixels (default 1.0; 0/None disables).
    symmetric_clim : montage display range. False (default) = ImageJ-style:
        lowest/highest saturation/2 percent of montage pixels clipped.
        True = [-m, +m] with m = (100 - saturation)th percentile of |pixels|.
    saturation : percent of montage pixels saturated (default 0.5); one global
        range for all tiles. vmin / vmax override it.
    mip_vmin, mip_vmax, mip_saturation : display range for the MIP panel. By
        default the floor is 0 (positive range only) and the ceiling is the
        (100 - mip_saturation)th percentile of MIP pixels (default 1.0%), which
        stretches contrast in the positive range.
    save_prefix : if given, writes
        <prefix>_montage.png         figure: montage + MIP (300 dpi)
        <prefix>_montage.tif         float32 stitched montage (true values)
        <prefix>_delta_images.tif    float32 stack, one image per event (N, H, W)
        <prefix>_location_images.tif float32 stack, one image per location,
                                     row by row (rows*cols, H, W)
        <prefix>_mip.tif             float32 maximum projection across locations
      The TIFFs hold the unscaled difference values (including negatives), not
      the display-scaled ones; ImageJ opens them as 32-bit images.
 
    Returns dict: delta_images, loc_images, montage, mip, xy, x_vals, y_vals,
    clim, mip_clim, fig, ax, ax_mip.
    """
    if stat not in ("mean", "max"):
        raise ValueError('stat must be "mean" or "max"')
    stack = np.moveaxis(np.asarray(stack), time_axis, 0)
    T, H, W = stack.shape
    events = np.atleast_2d(np.asarray(events, dtype=float))
    if events.shape[1] != 3:
        raise ValueError("events must have 3 columns: (intensity, start, stop)")
 
    n_ev, n_stim = len(events), len(stim_xy)
    starts = events[:, 1].astype(int)
    stops = events[:, 2].astype(int)
    win_ends = stops + post_frames  # last frame of each stimulus window
 
    delta = np.full((n_ev, H, W), np.nan, dtype=np.float32)
    for i, (s, e) in enumerate(zip(starts, stops)):
        s = max(s, 0)
        e_end = min(e + (1 if stop_inclusive else 0) + post_frames, T)
        b0 = max(s - baseline_frames, 0)
        if e_end <= s or s <= b0:
            print(f"Warning: event {i + 1} skipped (empty stimulus or baseline window).")
            continue
        if s - b0 < baseline_frames:
            print(f"Warning: event {i + 1} has only {s - b0} baseline frames available.")
        if np.any((starts < s) & (win_ends >= b0)):
            print(f"Warning: event {i + 1} baseline window overlaps an earlier "
                  f"stimulus window (incl. post_frames).")
 
        if stat == "mean":
            stim_img = stack[s:e_end].mean(axis=0, dtype=np.float32)
            if blur_sigma:
                stim_img = gaussian_filter(stim_img, blur_sigma)
        else:  # "max": blur each frame, then max over time
            w = np.asarray(stack[s:e_end], dtype=np.float32)
            if blur_sigma:
                w = gaussian_filter(w, sigma=(0, blur_sigma, blur_sigma))
            stim_img = w.max(axis=0)
        base_img = stack[b0:s].mean(axis=0, dtype=np.float32)
        if blur_sigma:
            base_img = gaussian_filter(base_img, blur_sigma)
        delta[i] = stim_img - base_img
 
    xy = np.array([stim_xy[i % n_stim] for i in range(n_ev)], dtype=float)
    x_vals = np.unique([p[0] for p in stim_xy])
    y_vals = np.unique([p[1] for p in stim_xy])[::-1]
    n_r, n_c = len(y_vals), len(x_vals)
 
    sums = np.zeros((n_r, n_c, H, W), dtype=np.float64)
    counts = np.zeros((n_r, n_c))
    for (x, y), img in zip(xy, delta):
        if np.isnan(img).all():
            continue
        xi = np.where(x_vals == x)[0][0]
        yi = np.where(y_vals == y)[0][0]
        sums[yi, xi] += img
        counts[yi, xi] += 1
    loc_images = np.full((n_r, n_c, H, W), np.nan, dtype=np.float32)
    for yi in range(n_r):
        for xi in range(n_c):
            if counts[yi, xi] > 0:
                loc_images[yi, xi] = sums[yi, xi] / counts[yi, xi]
 
    # (rows, cols, H, W) -> (rows*H, cols*W)
    montage = loc_images.transpose(0, 2, 1, 3).reshape(n_r * H, n_c * W)
 
    # Maximum intensity projection across all locations
    with np.errstate(all="ignore"):
        mip = np.nanmax(loc_images.reshape(-1, H, W), axis=0) \
            if np.isfinite(loc_images).any() else np.full((H, W), np.nan, np.float32)
 
    # Global display range for the montage (same LUT for every tile)
    finite = montage[np.isfinite(montage)]
    if vmin is None and vmax is None:
        if finite.size == 0:
            vmin, vmax = -1.0, 1.0
        elif symmetric_clim:
            m = np.percentile(np.abs(finite), 100 - saturation)
            vmin, vmax = -m, m
        else:
            vmin, vmax = np.percentile(finite, [saturation / 2, 100 - saturation / 2])
    elif vmin is None:
        vmin = float(finite.min()) if finite.size else -1.0
    elif vmax is None:
        vmax = float(finite.max()) if finite.size else 1.0
    if vmax <= vmin:
        vmax = vmin + 1e-12
 
    # Display range for the MIP (positive-range stretch)
    mfin = mip[np.isfinite(mip)]
    mv0 = 0.0 if mip_vmin is None else float(mip_vmin)
    if mip_vmax is None:
        mv1 = float(np.percentile(mfin, 100 - mip_saturation)) if mfin.size else mv0 + 1.0
    else:
        mv1 = float(mip_vmax)
    if mv1 <= mv0:
        mv1 = mv0 + 1e-12
 
    fig = ax = ax_mip = None
    if plot or save_prefix:
        top_h = 3 * n_c * (n_r * H) / (n_c * W)   # montage height (in) at 3 in/tile wide
        mip_h = 3.2
        fig = plt.figure(figsize=(1.8 + 3 * n_c, top_h + mip_h + 2.2),
                         layout="constrained")
        gs = fig.add_gridspec(2, 1, height_ratios=[top_h + 0.9, mip_h + 0.6])
        ax = fig.add_subplot(gs[0])
        ax_mip = fig.add_subplot(gs[1])
 
        im = ax.imshow(montage, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_xticks([(c + 0.5) * W - 0.5 for c in range(n_c)])
        ax.set_xticklabels([f"{v:g}" for v in x_vals])
        ax.set_yticks([(r + 0.5) * H - 0.5 for r in range(n_r)])
        ax.set_yticklabels([f"{v:g}" for v in y_vals])
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        for c in range(1, n_c):
            ax.axvline(c * W - 0.5, color="w", lw=0.5)
        for r in range(1, n_r):
            ax.axhline(r * H - 0.5, color="w", lw=0.5)
        fig.colorbar(im, ax=ax, label=f"Stimulus {stat} - baseline mean", shrink=0.8)
        ax.set_title(f"Response montage ({stat}; stim window +{post_frames} fr; "
                     f"baseline {baseline_frames} fr)")
 
        im2 = ax_mip.imshow(mip, cmap=cmap, vmin=mv0, vmax=mv1, interpolation="nearest")
        ax_mip.set_xticks([])
        ax_mip.set_yticks([])
        fig.colorbar(im2, ax=ax_mip, label="Max projection across locations", shrink=0.8)
        ax_mip.set_title("Maximum intensity projection (all locations)")
        ax.set_anchor("C")
        ax_mip.set_anchor("C")
 
    if save_prefix:
        d = os.path.dirname(save_prefix)
        if d:
            os.makedirs(d, exist_ok=True)
        fig.savefig(f"{save_prefix}_montage.png", dpi=300)
        import tifffile
        desc = (f"stat={stat}; baseline_frames={baseline_frames}; post_frames={post_frames}; "
                f"blur_sigma={blur_sigma}; units = stimulus {stat} minus baseline mean (raw values)")
        if Savetiff:
            tifffile.imwrite(f"{save_prefix}_montage.tif", montage.astype(np.float32), description=desc)
            tifffile.imwrite(f"{save_prefix}_delta_images.tif", delta.astype(np.float32), description=desc)
            tifffile.imwrite(f"{save_prefix}_location_images.tif",
                            loc_images.reshape(-1, H, W).astype(np.float32), description=desc)
            tifffile.imwrite(f"{save_prefix}_mip.tif", mip.astype(np.float32), description=desc)
        if not plot:
            plt.close(fig)
            fig = ax = ax_mip = None
 
    return {"delta_images": delta, "loc_images": loc_images, "montage": montage,
            "mip": mip, "xy": xy, "x_vals": x_vals, "y_vals": y_vals,
            "clim": (float(vmin), float(vmax)), "mip_clim": (mv0, mv1),
            "fig": fig, "ax": ax, "ax_mip": ax_mip}
 
 # ----------------------------------------------------------------------------

def otsu_threshold(values, nbins=256):
    """Otsu's threshold (maximises between-class variance) of a 1D set of values."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    hist, edges = np.histogram(v, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centers)
    with np.errstate(divide="ignore", invalid="ignore"):
        between = w0 * w1 * (m0 / w0 - (m0[-1] - m0) / w1) ** 2
    return centers[np.nanargmax(between)]
 
def make_roi_mask(mip, method="percentile", top_percent=10, cell_radius=5,
                  blob_k=5.0, min_area=15):
    """
    Turn a maximum-intensity-projection image into a labeled ROI mask.
 
    method
      "percentile": brightest `top_percent` percent of pixels (always selects
                    that fraction, whether or not there are cells).
      "otsu":       Otsu threshold on the MIP pixel histogram.
      "blob":       round-object detector. The MIP is filtered with a
                    Laplacian-of-Gaussian matched to `cell_radius` (sigma =
                    radius / sqrt(2)); local maxima of the response above
                    median + blob_k * (robust SD) are taken as cell centers, and
                    each cell is the connected region around its center that is
                    above half-maximum (peak vs. median background), limited to
                    1.5 x cell_radius from the center. Isolated bright pixels
                    and speckle are not round objects and are rejected.
    Regions smaller than `min_area` pixels are dropped.
 
    Returns (labels, info): labels is an int32 (H, W) image, 0 = background and
    1..K = one label per region/cell; info has the threshold and cell centers.
    """
    mip = np.asarray(mip, dtype=float)
    fin = np.isfinite(mip)
    if not fin.any():
        raise ValueError("MIP has no finite pixels")
    filled = np.where(fin, mip, np.median(mip[fin]))
    H, W = filled.shape
    info = {"method": method}
 
    if method in ("percentile", "otsu"):
        if method == "percentile":
            thr = np.percentile(filled[fin], 100 - top_percent)
            binary = filled >= thr
        else:
            thr = otsu_threshold(filled[fin])
            binary = filled > thr
        info["threshold"] = float(thr)
        lab, n = ndi.label(binary)
        if n:
            areas = ndi.sum(binary, lab, index=np.arange(1, n + 1))
            keep = np.concatenate([[False], areas >= min_area])
            lab, n = ndi.label(keep[lab])
        labels = lab.astype(np.int32)
    elif method == "blob":
        sigma = cell_radius / np.sqrt(2)
        resp = -sigma ** 2 * ndi.gaussian_laplace(filled, sigma)
        med = np.median(resp)
        mad = 1.4826 * np.median(np.abs(resp - med))
        thr = med + blob_k * mad
        size = int(2 * round(cell_radius) + 1)
        is_peak = (resp == ndi.maximum_filter(resp, size=size)) & (resp > thr)
        py, px = np.nonzero(is_peak)
        order = np.argsort(resp[py, px])[::-1]
        bg = np.median(filled)
        yy, xx = np.mgrid[:H, :W]
        labels = np.zeros((H, W), dtype=np.int32)
        n = 0
        for idx in order:
            y, x = py[idx], px[idx]
            pv = filled[y, x]
            if pv <= bg:
                continue
            half = bg + 0.5 * (pv - bg)
            cand = ((filled >= half) & (labels == 0) &
                    ((yy - y) ** 2 + (xx - x) ** 2 <= (1.5 * cell_radius) ** 2))
            cl, _ = ndi.label(cand)
            if cl[y, x] == 0:
                continue
            region = cl == cl[y, x]
            if region.sum() < min_area:
                continue
            n += 1
            labels[region] = n
        info["threshold"] = float(thr)
        info["n_peaks_found"] = int(len(py))
    else:
        raise ValueError('method must be "percentile", "otsu" or "blob"')
 
    n = int(labels.max())
    info["n_cells"] = n
    info["centers"] = (np.array(ndi.center_of_mass(labels > 0, labels, range(1, n + 1)))
                       if n else np.empty((0, 2)))  # (row, col)
    info["areas"] = (np.array(ndi.sum(labels > 0, labels, range(1, n + 1)))
                     if n else np.empty(0))
    return labels, info
 
def roi_trace(stack, mask, time_axis=0):
    """Mean of the raw (unblurred) stack over the pixels in `mask`, per frame."""
    stack = np.moveaxis(np.asarray(stack), time_axis, 0)
    return stack[:, np.asarray(mask, dtype=bool)].mean(axis=1, dtype=np.float64)
 
def analyze_masked_trace(stack, events, stim_xy=STIM_XY, mask_method="percentile",
                         top_percent=10, cell_radius=5, blob_k=5.0, min_area=15,
                         per_cell=False, baseline_frames=30, post_frames=10,
                         stop_inclusive=True, time_axis=0, stack_kwargs=None,
                         events_kwargs=None, plot=True, save_prefix=None):
    """
    1. analyze_stack  -> MIP across all stimulus locations.
    2. make_roi_mask  -> mask of the active area (see make_roi_mask for methods).
    3. roi_trace      -> mean raw trace inside the mask (T,).
    4. analyze_events -> peak / baseline-median ratio per stimulus + heat map.
 
    The same baseline_frames / post_frames are used for the image and trace
    analyses so the montage and the heat map are directly comparable.
    per_cell=True also runs step 3-4 separately for every labeled cell/region.
    stack_kwargs / events_kwargs: extra keyword arguments for analyze_stack /
    analyze_events (e.g. stat="max", blur_sigma=2, cmap=...).
 
    With save_prefix, in addition to the analyze_stack files, writes:
        <prefix>_mask.png         MIP with the mask outlines and cell numbers
        <prefix>_mask.tif         uint16 label image (0 = background)
        <prefix>_roi_trace.csv    mean trace of the whole mask
        <prefix>_roi_heatmap.png / _roi_raster.png / _roi_values.csv / _roi_grid.csv
                                  (analyze_events outputs for the pooled mask)
        <prefix>_cells_traces.csv, <prefix>_cells_ratio.csv   (per_cell only)
 
    Returns dict: stack (analyze_stack result), mip, labels, mask_info, trace,
    result (analyze_events result for the pooled mask), cells (per-cell
    analyze_events results keyed by label, or None), fig_mask.
    """
    stk = analyze_stack(stack, events, stim_xy=stim_xy, baseline_frames=baseline_frames,
                        post_frames=post_frames, stop_inclusive=stop_inclusive, blur_sigma=2,
                        time_axis=time_axis, plot=plot, save_prefix=save_prefix,
                        **(stack_kwargs or {}))
    mip = stk["mip"]
    labels, info = make_roi_mask(mip, method=mask_method, top_percent=top_percent,
                                 cell_radius=cell_radius, blob_k=blob_k,
                                 min_area=min_area)
    print(f"Mask ({mask_method}): {info['n_cells']} region(s), "
          f"{int((labels > 0).sum())} px")
    if info["n_cells"] == 0:
        print("Warning: no active regions found; no trace analysis was run.")
        return {"stack": stk, "mip": mip, "labels": labels, "mask_info": info,
                "trace": None, "result": None, "cells": None, "fig_mask": None}
 
    mask = labels > 0
    trace = roi_trace(stack, mask, time_axis)
    ev_kw = dict(baseline_frames=baseline_frames, post_frames=post_frames,
                 stop_inclusive=stop_inclusive)
    ev_kw.update(events_kwargs or {})
    result = analyze_events(trace, events, stim_xy=stim_xy, plot=plot,
                            save_prefix=f"{save_prefix}_roi" if save_prefix else None,
                            **ev_kw)
 
    cells = None
    if per_cell:
        cells, cell_traces = {}, []
        for lab in range(1, info["n_cells"] + 1):
            tr = roi_trace(stack, labels == lab, time_axis)
            cells[lab] = analyze_events(tr, events, stim_xy=stim_xy, plot=False,
                                        **ev_kw)
            cells[lab]["trace"] = tr
            cell_traces.append(tr)
 
    fig = None
    if plot or save_prefix:
        fig, ax = plt.subplots(figsize=(6, 6 * mip.shape[0] / mip.shape[1] + 0.6))
        im = ax.imshow(mip, cmap="magma", vmin=stk["mip_clim"][0], vmax=stk["mip_clim"][1])
        ax.contour(mask.astype(float), levels=[0.5], colors="cyan", linewidths=1)
        for lab, (cy, cx) in enumerate(info["centers"], start=1):
            ax.text(cx, cy, str(lab), color="w", ha="center", va="center", fontsize=1)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"ROI mask ({mask_method}): {info['n_cells']} region(s), {int(mask.sum())} px")
        fig.colorbar(im, ax=ax, shrink=0.8, label="Max projection across locations")
        fig.tight_layout()
 
    if save_prefix:
        import tifffile
        d = os.path.dirname(save_prefix)
        if d:
            os.makedirs(d, exist_ok=True)
        fig.savefig(f"{save_prefix}_mask.png", dpi=300)
        #tifffile.imwrite(f"{save_prefix}_mask.tif", labels.astype(np.uint16))
        np.savetxt(f"{save_prefix}_FOV_trace.csv", trace, delimiter=",", fmt="%.6g",
                   header="mean_trace", comments="")
        if per_cell:
            cols = [f"cell{k}" for k in cells]
            np.savetxt(f"{save_prefix}_cells_traces.csv", np.column_stack(cell_traces),
                       delimiter=",", fmt="%.6g", header=",".join(cols), comments="")
            #np.savetxt(f"{save_prefix}_cells_ratio.csv",
            #           np.column_stack([cells[k]["ratio"] for k in cells]),
            #           delimiter=",", fmt="%.6g", header=",".join(cols), comments="")
        if not plot:
            plt.close(fig)
            fig = None
 
    return {"stack": stk, "mip": mip, "labels": labels, "mask_info": info,
            "trace": trace, "result": result, "cells": cells, "fig_mask": fig}

##Full Main program
Meta_Data = extract_TS_metadata() ##extract metadata of image files
F,MaxZTiff = Load_OMEtiff(Meta_Data) #use metadata to load the OME Tiffs and rearrange the images and save a MaxZ for use in analysis
Frame_Force = Load_ThorSync(flybackCycles=Meta_Data["flybackCycles"]) ##Load thorSync file for finding stimulus peaks. flyback is (Z+flyback) as frame counter includes flybacks
mN_Start_Stop = Find_Peaks(Frame_Force, Imaging_Folder) ##find the frame start and stops of the mech stim
analyze_stack(MaxZTiff, mN_Start_Stop, save_prefix = f"{Imaging_Folder}RF",Savetiff=False) #Make image montage of each stimulus and arrange into RF map.
out = analyze_masked_trace(MaxZTiff, mN_Start_Stop, mask_method="blob", top_percent=10, cell_radius=5, blob_k = 5, per_cell=True, save_prefix=f"{Imaging_Folder}RF")
##Create masks of neurons and extract fluoresence traces to then find the peak responses to each stimulus and map into a RF map.
##Note on methods for masks: 
    ## Decrease the k to increase the number of blobs considered as ROIS. Cell_radius also strongly gates what is selected
    ##Alternatively use percentile which doesn't care what shape the response comes in