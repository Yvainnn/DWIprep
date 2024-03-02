 #!/bin/bash

# ____________________DWI preprocessing_____________________  

#-------------The folders will be organized as: ------------
#		
# PRJ_main folder 
#------	SBJ_n main folder		-
#	      sess01
#			anat
#				freesurfer_output
#			diff
#				01_denoise		
#				02_degibbs
#				03_topup
#				04_eddy
#				05_biascorrect
#				06_dtifit
#				07_diff2std
#			t2
#				01_motioncorrect
#				02_curvefitting
#				03_distorcorrect
#				04_biascorrect&Jacobian
#				05_qt2std
#		sess02
#			anat
#				freesurfer_output
#			diff
#				01_denoise		
#				02_degibbs
#				03_topup
#				04_eddy
#				05_biascorrect
#				06_dtifit
#				07_diff2std
#			t2
#				01_motioncorrect
#				02_curvefitting
#				03_distorcorrect
#				04_biascorrect&Jacobian
#				05_qt2std
#		sess03
#		...... (session folders)
#		
#		std2mni_matrix
#		sess012sess0n
#
#
#--------------------------------------------------------------


#		0) -- Bids setup and dicom2niix (better run it before start processing or check it step by step)

dicom_path=/home ecc...
diff_scan= #name of dicom diff scans
anat_scan= #name of dicom anat/t1 scans
dcm2niix_path=/home/malberti/Unix_Folders/dcm2niix
prj="IDK" #project name
session= #session number
subj= #name of working directory


ls "$dicom_path" > subj_list.txt #This line creates a txt file containing all the name of all the subj of the project

for subj in `cat subj_list.txt`;
do
	mkdir "$subj"
	cd "$subj"
	for i in {1..$session};		#This loop should set up all the folder that the dwi processing will needs
	do
		mkdir sess0"$session"
		cd sess0"$session"
		mkdir anat
		mkdir diff
		cd diff
		mkdir 01_denoise		
		mkdir 02_degibbs
 		mkdir 03_topup
		mkdir 04_eddy
		mkdir 05_biascorrect
		mkdir 06_dtifit
		mkdir 07_diff2std
		mkdir 08_csd
		mkdir 09_tractography
		
		cd /home/malberti/Unix_Folders/"$prf"
		
	done

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# (2) dicom2niix --> 

./dicom2niix_convert.sh

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# (3) DWI preprocessing

nsub= # number of subjects
subj_prefix= # subjects' common prefix for folders
nsession= # number of sessions

workig_dir= # main path; this folder should contain all the subject folders

for ((sbj=0; sbj<="$nsub"; sbj++)) # moving through subjects
do 
    for ((nss=1; nss<="$nsession"; nss++)) 
    do
        # each session contains 4 folders named as diff_30_AP/PA and diff_96_AP/PA. Denoise and Gibbs ring correction should be applied to each sequence. To do it, you can write a file named 'dwi_raw.txt':
        # ---------------------
        # diff_30_PA
        # diff_30_PA
        # diff_96_AP
        # diff_96_PA
        # ------------------\
        sess_path="$working_dir"/"$subj_prefix"_"$sbj"/sess0"$nss" # this folder should contain all the folders created above and the raw NIfTI file
        denoise="$sess_path"/01_denoise
        gibbs="$sess_path"/02_degibbs
        topup="$sess_path"/03_topup
        eddy="$sess_path"/04_eddy
        biascorrect="$sess_path"/05_biascorrect
        dtifit="$sess_path"/06_dtifit

        acqparams= # add path of acq params file 	
        index= # add path of index file 
		for nifti in `cat dwi_raw.txt`
		do 
			dwi="$sess_path"/"$nifti"
           		mrconvert -fslgrad "$dwi"/diff.bvec "$dwi"/diff.bval "$dwi"/dwi.nii.gz  "$dwi"/dwi.mif # convert NIfTI file to MRtrix3 MIF format
            		dwiextract "$dwi"/dwi.mif "$dwi"/dwi_b0.mif -bzero -force # extract all b0 images
           		mrconvert "$dwi"/dwi_b0.mif "$dwi"/dwi_b0.nii.gz 
            		fslmaths "$dwi"/dwi_b0.nii.gz -Tmean "$denoise"/"$nifti"_b0.nii.gz # create the mean b0
            		bet  "$denoise"/"$nifti"_b0.nii.gz  "$denoise"/"$nifti"_b0 -R -f 0.1 -g 0 -m # create a brain mask from the b0

            		# [Denoise]
            		dwidenoise -estimator Exp2 -mask "$denoise"/"$nifti"_b0_mask.nii.gz "$dwi"/dwi.mif "$denoise"/"$nifti"_denoised_diff.mif \
	      				-noise "$denoise"/"$nifti"_noise.mif -force # run denoising
            		mrcalc "$dwi"/dwi.mif "$denoise"/"$nifti"_denoised_diff.mif -subtract "$denoise"/"$nifti"_diff_res.nii.gz -force # estimate the difference between raw DWI and denoised

            		# [Gibbs ring correction]
            		mrdegibbs "$denoise"/"$nifti"_denoised_diff.mif  "$gibbs"/"$nifti"_denoised_gibbs.mif -force # run Gibbs ring correction
            		mrcalc "$denoise"/"$nifti"_denoised_diff.mif  "$gibbs"/"$nifti"_denoised_gibbs.mif -subtract "$gibbs"/"$nifti"_gibbs_res.mif # estimate the difference between denoised and denoised_degibbs images
            		mrconvert "$gibbs"/"$nifti"_denoised_gibbs.mif "$gibbs"/"$nifti"_denoised_gibbs.nii.gz

            		# Topup and eddy correction will be performed on DWI images obtained by merging 30_AP/PA and 96_AP/PA b0s. 
            		dwiextract "$gibbs"/"$nifti"_denoised_gibbs.mif "$topup"/b0_"$nifti".mif -bzero -force # get the corrected b0 shell 
            		mrconvert "$topup"/b0_"$nifti".mif "$topup"/b0_"$nifti".nii.gz
            		paste "$nifti"/dwi.bvec >> "$eddy"/allbvecs.bvec # create all bvecs and all bval files 
            		paste "$nifti"/dwi.bval >> "$eddy"/allbval.bval
        	done

        	# Now go back to subj/sess"n" folder, merge everything and run topup -> eddy -> bias correct

        	# [Topup]		
        	fslmerge -t "$eddy"/alldwi.nii.gz "$gibbs"/diff*_denoised_gibbs.nii.gz # merge all DWI and all b0 images  
        	fslmerge -t "$topup"/b0_images.nii.gz "$topup"/b0_diff*.nii.gz 
        	# CHECK THE COHERENCE BETWEEN DWI<->BVEC<->BVAL and if alldwi.nii.gz is coherent with acqparams/index
        	topup --imain="$topup"/b0_images.nii.gz \
            		--datain="$acqparams" --config=/usr/local/fsl/etc/flirtsch/b02b0.cnf \
            		--out="$topup"/topup --fout="$topup"/topup_rfields \
            		--iout="$topup"/topup_unwarped --verbose # run topup on merged b0

        	fslmaths "$topup"/topup_unwarped.nii.gz -Tmean "$eddy"/unwarped_b0_mean.nii.gz # create a new brain mask from topup unwarped output
        	bet "$eddy"/unwarped_b0_mean.nii.gz "$eddy"/Tmean_b0 -R -f 0.1 -g 0 -m

        	# [Eddy current]
        	eddy --imain="$eddy"/alldwi.nii.gz --mask="$eddy"/Tmean_b0_mask.nii.gz \
            		--acqp="$acqparams" --index="$index" --bvecs="$eddy"/allbvecs.bvec \
            		--bvals="$eddy"/allbval.bval --topup="$topup"/topup \
            		--out="$eddy"/eddy_out --flm=quadratic --slm=linear \
            		--dont_peas --verbose # run eddy current

        	# [Bias correction] 
        	dwibiascorrect -ants "$eddy"/eddy_out.nii.gz \
            			"$biascorrect"/eddy_unbiased_subj"$sbj"_sess0"$nss".mif \
            			-fslgrad "$eddy"/eddy_out.eddy_rotated_bvecs "$eddy"/allbval.bval \
            			-bias "$biascorrect"/subj"$sbj"_sess0"$nss"_bias.mif # run bias correction

        	mrconvert "$biascorrect"/eddy_unbiased_subj"$sbj"_sess0"nss".mif "$biascorrect"/eddy_unbiased_subj"$sbj"_sess0"nss".nii.gz # convert processed DWI from MIF to NIfTI

        	echo "sub: "$sbj", session: 0"$nss" Done successfully"
    	done
done	
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Next steps should be modified based on what you want to compute or the type of DWI sequence you are processing.
# If necessary, you can now run DTIfit to estimate the tensor, FA, etc. Diffusivity metrics can also be estimated using fODF later in the pipeline.
# Add these lines to the previous script:
	dwiextract "$biascorrect"/eddy_unbiased_subj"$sbj"_sess0"nss".mif "$dtifit"/subj"$sbj"_sess0"nss"_b1000.mif  \
				-export_grad_fsl  "$dtifit"/b1000_bvec.bvec  "$dtifit"/b1000_bval.bval  -singleshell 1000 -force # extract the b1000 shell to compute the tensor; check if it works because sometimes it fails to extract b995 or b1005 shells
	mrconvert "$dtifit"/subj"$sbj"_sess0"nss"_b1000.mif  "$dtifit"/subj"$sbj"_sess0"nss"_b1000.nii.gz
	dtifit -k "$dtifit"/subj"$sbj"_sess0"nss"_b1000.nii.gz -m "$eddy"/Tmean_b0_mask.nii.gz -o "$dtifit"/subj"$sbj"_sess0"nss" \
			-r "$dtifit"/b1000_bvec.bvec -b "$dtifit"/b1000_bval.bval -V --save_tensor  # run FSL dtifit on b1000 shell

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# DTIFIT - if you want to compute multi-shell/kurtosis data (b0, b1000, and b2500) [I haven't tested it yet]
	mkdir "$sess_path"/08_kurtosis
	kurtosis="$sess_path"/08_kurtosis
	dwiextract "$biascorrect"/eddy_unbiased_subj"$sbj"_sess0"nss".mif "$kurtosis"/subj"$sbj"_sess0"nss"_kurtosis.mif \
				-export_grad_fsl  "$kurtosis"/kurtosis_bvec.bvec  "$kurtosis"/kurtosis_bval.bval  -singleshell 1000,2500 -force  # extract b0, b1000, and b2500 shells [Check it]
	mrconvert "$kurtosis"/subj"$sbj"_sess0"nss"_kurtosis.mif  "$kurtosis"/subj"$sbj"_sess0"nss"_kurtosis.nii.gz
	dtifit -k "$kurtosis"/subj"$sbj"_sess0"nss"_kurtosis.nii.gz  -m "$eddy"/Tmean_b0_mask.nii.gz -o "$kurtosis"/k_subj"$sbj"_sess0"nss" \
			-r "$kurtosis"/kurtosis_bvec.bvec -b "$kurtosis"/kurtosis_bval.bval -V --save_tensor  --kurt --kurtdir   # run FSL dtifit on multi-shell data

# -------------------------------------------------------------------------------------------------------------------------------------" ---------------------------------------------------

# Run TBSS analysis. To run it, you need the DWI metrics map obtained from DTIfit.
# Note: This script runs outside the main loop. It may also work with DKI images, but it needs to be tested.
# (1) - TBSS preproc should work even without copying all metrics files into the same folder. 
# Usually, I copy everything into the same folder, so...
# Anyway, set "$working_dir" before running.

mkdir "$working_dir"/tbss
tbss="$working_dir"/tbss

mkdir "$tbss"/FA #directory and folders setup
#mkdir "$tbss"/MD
#mkdir "$tbss"/AD
#mkdir "$tbss"/RD

fa="$tbss"/FA
#md="$tbss"/MD #To run it with MD just change the input file, dtifit compute it by itself
#ad="$tbss"/AD #To compute AD tbss run it on L1 files, dtifit compute it 
#rd="$tbss"/RD #To run RD tbss -> RD has to be computed from dtifit output, the RD is the average between L2 and L3 files so: 

# compute RD
# fslmerge -T temp_mean.nii.gz L3.nii.gz L2.nii.gz #do it for each subject
# fslmaths temp_mean.nii.gz -Tmean RD.nii.gz #compute the mean between L2 and L3 and obtain the radial diffusivity map

# to set up the right order u have to create one array that contain the bname of every images in the right order to do it just:

dwi_name=()
for ((nss=1;nss<="$nsession";nss++)) 
do 
	for ((sbj=0;sbj<="$nsub";sbj++)) #moving throught subject
	do 
		sess_path= "$working_dir"/"$subj_prefix"_"$subj"/sess0"$nss"
 		dtifit="$sess_path"/06_dtifit
		dwi_name+=""$dtifit"/subj"$sbj"_sess0"nss"_FA.nii.gz"
	done

done 

echo "${dwi_name}"
cd "$fa" #check the order before running tbss, but it should be

tbss_1_preproc "${dwi_name}" 
tbss_2_reg  -T
#Target-selection options - choose ONE of:
# -T            : use FMRIB58_FA_1mm as target for nonlinear registrations (recommended)
# -t <target>   : use <target> image as target for nonlinear registrations
# -n            : find best target from all images in FA

tbss_3_postreg -S
#Choose ONE of:
# -S   : derive mean_FA and mean_FA_skeleton from mean of all subjects in study (recommended)
# -T   : use FMRIB58_FA and its skeleton instead of study-derived mean and skeleton

tbss_4_prestats 0.2 

#Now the design matrices have to be created, the run randomise 

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## Next steps have been elaborated from DWI course UNIL-FENS 2024
# Tractography, structural connectivity, and fODF estimation + computing tensor from fODF [multishell DWI data] - [it hasn't been tested yet]

for ((sbj=0;sbj<="$nsub";sbj++)) # moving through subject
do 
	for ((nss=1;nss<="$nsession";nss++)) 
	do
		sess_path="$working_dir"/"$subj_prefix"_"$subj"/sess0"$nss" # this folder should contain all the folders created above and the raw NII file
		biascorrect="$sess_path"/05_biascorrect
		csd="$sess_path"/08_csd
		anat="$sess_path"/anat/ # add the name of processed T1 images
		tract="$sess_path"/09_tractography 
		n_tract= # number of streamlines that you want to compute in tckgen
		reduced_tract= # number of tracts that you want to obtain by tcksift

		bet "$biascorrect"/eddy_unbiased_subj"$sbj"_sess0"nss".nii.gz "$csd"/subj"$sbj"_sess0"nss"_brain -R -f 0.2 -m -g 0 # create a mask from the last pre-processing output [dwi2mask could be a good option too]
		mrconvert "$csd"/subj"$sbj"_sess0"nss"_brain_mask.nii.gz "$csd"/subj"$sbj"_sess0"nss"_brain_mask.mif
		mrconvert "$csd"/subj"$sbj"_sess0"nss"_brain.nii.gz "$csd"/subj"$sbj"_sess0"nss"_brain.mif
		dwiextract "$csd"/subj"$sbj"_sess0"nss"_brain.mif "$csd"/subj"$sbj"_sess0"nss"_b0.mif -bzero -force # extract corrected b0 and create a mean b0 image
		mrmath "$csd"/subj"$sbj"_sess0"nss"_b0.mif mean -axis 3 "$csd"/subj"$sbj"_sess0"nss"_b0.mif -force

#		[Optional] Run interpolation
#
#		mrgrid "$csd"/subj"$sbj"_sess0"nss"_brain.mif regrid "$csd"/subj"$sbj"_sess0"nss"_brain_reg.mif -voxel 1.5 # Interpolating the diffusion data also interpolates the data to be used on the diffusion images (anatomical plots and brain mask).
#		mrgrid "$csd"/subj"$sbj"_sess0"nss"_brain_mask.mif regrid "$csd"/subj"$sbj"_sess0"nss"_brain_mask_reg.mif -voxel 1.5 -interp nearest
#
#
		dwi2response dhollander "$csd"/subj"$sbj"_sess0"nss"_brain.mif "$csd"/response_wm.txt "$csd"/response_gm.txt "$csd"/response_csf.txt \
				  -mask "$csd"/subj"$sbj"_sess0"nss"_brain_mask.mif \
				  -voxels "$csd"/voxel.mif -force # compute response function
#				  -lmax LMAX # The maximum harmonic degree(s) for response function estimation (comma-separated list in case of multiple b-values)
#		shview <response> # to check the response function

# 	Estimation using CSD for multi-shell data; this can be achieved using the MSMT-CSD (Multi-Shell Multi Tissue CSD)
		dwi2fod msmt_csd "$csd"/subj"$sbj"_sess0"nss"_brain.mif \
				"$csd"/response_wm.txt "$csd"/wmfod.mif \
				"$csd"/response_gm.txt "$csd"/gm.mif \
				"$csd"/response_csf.txt "$csd"/csf.mif \
				-mask "$csd"/dwi_mask_upsampled.mif -force # estimate FOD for each tissue
#   			-lmax order
# 	The maximum spherical harmonic order for the output FOD(s). For algorithms with multiple outputs, this should be provided as a comma-separated list of integers, 
#	one for each output image; for single-output algorithms, only a single integer should be provided. 
#	If omitted, the command will use the lmax of the corresponding response function (i.e., based on its number of coefficients), up to a maximum of 8.

		mrconvert "$csd"/wmfod.mif "$csd"/wm.mif -coord 3 0 # extract the volume fraction
		mrcat "$csd"/csf.mif "$csd"/gm.mif "$csd"/wm.mif "$csd"/volume_fraction.mif # concatenate each volume fraction to check it

		mtnormalise "$csd"/wmfod.mif "$csd"/wmfod.mif \
					"$csd"/gm.mif "$csd"/gm.mif \
					"$csd"/csf.mif csf.mif -mask \
					"$csd"/dwi_mask_upsampled.mif -force # normalize intensity [OPTIONAL]

#	In the next line, check if the anat image is as NIFTI or MIF 
		5ttgen fsl -sgm_amyg_hipp "$anat" "$tract"/5tt_nocoreg.mif # generate a five-tissue-type (5TT) segmented tissue image suitable for use in Anatomically-Constrained Tractography (ACT)

# Intra-subj/session registration b02t1w
		mrconvert "$csd"/subj"$sbj"_sess0"nss"_b0.mif "$csd"/subj"$sbj"_sess0"nss"_b0.nii.gz
		flirt -in "$csd"/subj"$sbj"_sess0"nss"_b0.nii.gz -ref "$anat" -dof 6 -omat "$tract"/b02t1_fsl.mat # this line creates an FSL transformation matrix; you should convert it as an MRtrix3 matrix (next line)
		transformconvert "$tract"/b02t1_fsl.mat "$csd"/subj"$sbj"_sess0"nss"_b0.nii.gz "$anat" flirt_import "$tract"/b02t1_mrtrix.txt -force 

#	Intra-subj/session registration t1w2dwi 
		mrtransform "$anat" -linear "$tract"/b02t1_mrtrix.txt -inverse "$tract"/T1_coreg.mif -force

#	Intra-subj/session registration 5tt2dwi
		mrtransform "$tract"/5tt_nocoreg.mif -linear "$tract"/b02t1_mrtrix.txt -inverse "$tract"/5tt_coreg.mif -force
		
# Tractography		
		5tt2gmwmi "$tract"/5tt_coreg.mif "$tract"/gmwmSeed.mif -force
		tckgen -act "$tract"/5tt_coreg.mif -backtrack -crop_at_gmwmi -seed_gmwmi "$tract"/gmwmSeed.mif -select "$n_tract" "$csd"/wmfod.mif "$tract"/subj"$sbj"_sess0"nss".tck # run tractography reconstruction
		tcksift -act "$tract"/5tt_coreg.mif -remove_untracked -term_number "$reduced_tract" "$tract"/subj"$sbj"_sess0"nss".tck "$csd"/wmfod.mif "$tract"/subj"$sbj"_sess0"nss"_red.tck
		
		tck2connectome -symmetric -zero_diagonal -scale_invnodevol "$tract"/subj"$sbj"_sess0"nss".tck "$anat"/hcpmmp1_parcels_coreg.mif "$tract"/matrix_subj"$sbj"_sess0"nss".csv -out_assignment "$tract"/assignments_subj"$sbj"_sess0"nss".csv # "unweighted" connectivity matrix'

# DWI metrics weighted connectivity matrix
# To run it, you need the maps of the metric that you want to apply as weight to your structural connectivity matrix. You can use the maps created in dtifit lines and convert them into MIF format or estimate them again by using tensor2metric and dwi2tensor
# try to obtain the tensor metrics by using fod2dec, in order to get the FA map from fODF 

        dwi2tensor "$csd"/subj"$sbj"_sess0"nss"_brain.mif -mask "$csd"/subj"$sbj"_sess0"nss"_brain_mask.mif "$tract"/tensor.mif -force # estimate the tensor 
		tensor2metric "$tract"/tensor.mif -fa "$tract"/fa.mif -force # tensor2metric can estimate several diffusivity metrics; it doesn't work with multishell data 
		tcksample "$tract"/subj"$sbj"_sess0"nss".tck "$tract"/fa.mif "$tract"/subj"$sbj"_sess0"nss"_FA.csv -stat_tck mean -precise -use_tdi_fraction
		tck2connectome "$tract"/subj"$sbj"_sess0"nss".tck "$anat"/hcpmmp1_parcels_coreg.mif "$tract"/matrix_subj"$sbj"_sess0"nss"_FA.csv \
						-scale_file "$tract"/subj"$sbj"_sess0"nss"_FA.csv -stat_edge mean -symmetric -zero_diagonal -out_assignment "$tract"/assignments_subj"$sbj"_sess0"nss".csv


	done
done 


#	population_template	Generates an unbiased group-average template from a series of images
#	tckdfc	Perform the Track-Weighted Dynamic Functional Connectivity (TW-dFC) method
# 	 tckdfc [ options ] tracks fmri output
#
#        tracks       the input track file.
#        fmri         the pre-processed fMRI time series=
#        output       the output TW-dFC image

#     This command generates a Track-Weighted Image (TWI), where the  contribution from each streamline to the image is the Pearson correlation between the fMRI time series at the streamline endpoints.
#	  The output image can be generated in one of two ways (note that one of
#     these two command-line options MUST be provided): 

#     - "Static" functional connectivity (-static option): Each streamline
#     contributes to a static 3D output image based on the correlation between
#     the signals at the streamline endpoints using the entirety of the input
#     time series.

#     - "Dynamic" functional connectivity (-dynamic option): The output image is
#     a 4D image, with the same number of volumes as the input fMRI time series.
#     For each volume, the contribution from each streamline is calculated based
#     on a finite-width sliding time window, centred at the timepoint
#     corresponding to that volume.

#     Note that the -backtrack option in this command is similar, but not
#     precisely equivalent, to back-tracking as can be used with
#     Anatomically-Constrained Tractography (ACT) in the tckgen command.
#     However, here the feature does not change the streamlines trajectories in
#     any way; it simply enables detection of the fact that the input fMRI image
#     may not contain a valid timeseries underneath the streamline endpoint, and
#     where this occurs, searches from the streamline endpoint inwards along the
#     streamline trajectory in search of a valid timeseries to sample from the
#     input image.
