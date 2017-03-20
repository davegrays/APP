#!/bin/bash

reflist="DGr_NIMS_MRI_MSEL_ADOS_ADI_Query_20150610_caseIDonly.txt"
secondary_reflist="DGr_Connectome_BEH_20160420.txt"
MOI="DQ" #MOI is your Measure-of-Interest from the secondary_reflist
cfolder="../mapped/connectomes_diffusion_MRtrixCSD_prob"
topdir="../mapped"

#create list of finished subjects from the connectomes folder, and get number of subjects
echo "creating list of finished subjects from the connectomes folder"
numids=`ls $cfolder | wc -l`
ls $cfolder > NBSprep/SUBS_${numids}ids.txt

#scroll through subject list to create the corresponding design matrix based on subject diagnosis
#diagnoses must be set to ASD or Typical. This script is designed to exit with error if the subject list includes some of the subjects with other diagnoses (i.e. "See_Comment" or "ASD-FXS").
echo "creating corresponding design matrix."
echo "using subject information in ${reflist}."
echo "checking to make sure there are no subjects with diagnoses other than ASD or Typical"
rm -f NBSprep/design${numids}.txt
rm -f NBSprep/designANCOVA${numids}.txt
rm -f NBSprep/designANCOVA_withdiagXgender_${numids}.txt
rm -f NBSprep/designANCOVA_with${MOI}_${numids}.txt

for sub in `cat NBSprep/SUBS_${numids}ids.txt`;do
	fields=`cat ${reflist} | grep "$sub	"`
	diag=`echo $fields | awk '{print $3}'`
	age=`echo $fields | awk '{print $6}'`
	ICV=`cat ${topdir}/${sub}/FREESURFER/stats/aseg.stats  | grep ICV | awk '{print $7}' | sed 's/\,//'`
	gender=`echo $fields | awk '{print $2}'`

	wildcard=`echo $fields | awk '{print $2}'`
	if [ $gender -eq 1 ];then gender=1
	elif [ $gender -eq 2 ];then gender=-1
	else echo "Gender value $gender for subject $sub is invalid. Exiting early.";fi

	if [ "$diag" == "ASD" ];then
		echo "0 1" >> NBSprep/design${numids}.txt
		echo "1 -1 $gender $age $ICV" >> NBSprep/designANCOVA${numids}.txt
		echo "1 -1 $gender $((-1 * $gender)) $age $ICV" >> NBSprep/designANCOVA_withdiagXgender_${numids}.txt
		echo "1 -1 $gender $MOI $age $ICV" >> NBSprep/designANCOVA_with${MOI}_${numids}.txt
	elif [ "$diag" == "Typical" ];then
		echo "1 0" >> NBSprep/design${numids}.txt
		echo "1 1 $gender $age $ICV" >> NBSprep/designANCOVA${numids}.txt
		echo "1 1 $gender $gender $age $ICV" >> NBSprep/designANCOVA_withdiagXgender_${numids}.txt
		echo "1 1 $gender $MOI $age $ICV" >> NBSprep/designANCOVA_with${MOI}_${numids}.txt
	else
		echo "Diagnosis for sub ${sub} (given as ${diag}) is invalid. Must be ASD or Typical."
		echo "Removeing sub ${sub} from final list."
		sed -i "/${sub}/d" NBSprep/SUBS_${numids}ids.txt
		prevnum=${numids}
		numids=$(($numids-1))
		echo "Changing number of subjects from $prevnum to $numids"
		mv NBSprep/design${prevnum}.txt NBSprep/design${numids}.txt
		mv NBSprep/designANCOVA${prevnum}.txt NBSprep/designANCOVA${numids}.txt
		mv NBSprep/SUBS_${prevnum}ids.txt NBSprep/SUBS_${numids}ids.txt
		mv NBSprep/designANCOVA_withdiagXgender_${prevnum}.txt NBSprep/designANCOVA_withdiagXgender_${numids}.txt
		mv NBSprep/designANCOVA_with${MOI}_${prevnum}.txt NBSprep/designANCOVA_with${MOI}_${numids}.txt
	fi
done

echo "making group list for matlab matrices2multidim.m script"
echo "coded as TYP-m (1), TYP-f (2), ASD-m (3), ASD-f (4)"
cp -f NBSprep/designANCOVA_withdiagXgender_${numids}.txt explore/groups_${numids}ids.txt
sed -i 's/^1 1 1 1.*/1/' explore/groups_${numids}ids.txt
sed -i 's/^1 1 -1 -1.*/2/' explore/groups_${numids}ids.txt
sed -i 's/^1 -1 1 -1.*/3/' explore/groups_${numids}ids.txt
sed -i 's/^1 -1 -1 1.*/4/' explore/groups_${numids}ids.txt

echo "making list of genders for ASD group"
sed '/1/d' explore/groups_${numids}ids.txt | sed '/2/d' > explore/groups_${numids}ids_ASD_genders.txt
echo "making list of genders for TYP group"
sed '/3/d' explore/groups_${numids}ids.txt | sed '/4/d' > explore/groups_${numids}ids_TYP_genders.txt
echo "making list of diagnoses for males"
sed '/2/d' explore/groups_${numids}ids.txt | sed '/4/d' > explore/groups_${numids}ids_male_diagnoses.txt
echo "making list of diagnoses for females"
sed '/1/d' explore/groups_${numids}ids.txt | sed '/3/d' > explore/groups_${numids}ids_female_diagnoses.txt

echo "starting up the python script to do some actual number crunching..."
./printGroupInfo.py $numids

exit
