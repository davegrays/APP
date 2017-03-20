#!/usr/global/bin/python2.7

import os, errno, argparse

#recognize command line argument (# of subjects) 
parser = argparse.ArgumentParser()
parser.add_argument("number", help="the number of subjects listed in the design file name", type=int)
args = parser.parse_args()
number_of_subjects=args.number

#other I/O, hardcoded
fileprefix="NBSprep/designANCOVA"
filepreprefix="NBSprep/design"
subjectID_fileprefix="NBSprep/SUBS_"
outdir="explore"

#import number-crunching libraries
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def rm_f(filename):
    if os.path.exists(filename):
	try:
		os.remove(filename)
	except OSError, e:
		print ("Error: %s - %s." % (e.filename,e.strerror))

def get_data(fileprefix,number_of_subjects,subjectID_fileprefix):
    #load data
    data=np.loadtxt(fileprefix + str(number_of_subjects) + '.txt')
    subIDs=np.loadtxt(subjectID_fileprefix + str(number_of_subjects) + 'ids.txt',dtype=str)
    #define each column number (minus 1) as a metric
    metric = {
        'diagnosis' : 1,
        'gender' : 2,
        'age' : 3,
        'TCV' : 4,
    }
    return data, subIDs, metric

def make_gender_lists(fileprefix,filepreprefix,subjectID_fileprefix,number_of_subjects):
    #get data
    data, subIDs, metric = get_data(fileprefix,number_of_subjects,subjectID_fileprefix)
    diag_status=np.loadtxt(filepreprefix + str(number_of_subjects) + '.txt')
    #define male and female groups separately
    maleyes=data[:,metric['gender']]>0
    femaleyes=data[:,metric['gender']]<0
    #remove gender column
    data=np.delete(data,metric['gender'],1)
    #retrieve data
    maledata=data[maleyes]
    femaledata=data[femaleyes]
    male_diag_status=diag_status[maleyes]
    female_diag_status=diag_status[femaleyes]
    #retrieve IDs
    maleIDs=subIDs[maleyes]
    femaleIDs=subIDs[femaleyes]
    #print full ANCOVA data to txt files
    np.savetxt(fileprefix + str(number_of_subjects) + '_' + str(len(maledata)) + 'male.txt', maledata, fmt = ["%.0f"]*2 + ["%.1f"] + ["%.6f"])
    np.savetxt(fileprefix + str(number_of_subjects) + '_' + str(len(femaledata)) + 'female.txt', femaledata, fmt = ["%.0f"]*2 + ["%.1f"] + ["%.6f"])
    #print diagnosis only data to txt files
    np.savetxt(filepreprefix + str(number_of_subjects) + '_' + str(len(male_diag_status)) + 'male.txt', male_diag_status, fmt = ["%.0f"]*2)
    np.savetxt(filepreprefix + str(number_of_subjects) + '_' + str(len(female_diag_status)) + 'female.txt', female_diag_status, fmt = ["%.0f"]*2)
    #print IDs to txt files
    np.savetxt(subjectID_fileprefix + str(number_of_subjects) + '_' + str(len(maleIDs)) + 'maleids.txt', maleIDs, fmt = ["%s"])
    np.savetxt(subjectID_fileprefix + str(number_of_subjects) + '_' + str(len(femaleIDs)) + 'femaleids.txt', femaleIDs, fmt = ["%s"])

def make_diag_lists(fileprefix,subjectID_fileprefix,number_of_subjects):
    #get data
    data, subIDs, metric = get_data(fileprefix,number_of_subjects,subjectID_fileprefix)
    #define ASD and TYP groups separately
    TYPyes=data[:,metric['diagnosis']]>0
    ASDyes=data[:,metric['diagnosis']]<0
    #remove diagnosis column
    data=np.delete(data,metric['diagnosis'],1)
    #retrieve data
    TYPdata=data[TYPyes]
    ASDdata=data[ASDyes]
    #retrieve IDs
    TYPIDs=subIDs[TYPyes]
    ASDIDs=subIDs[ASDyes]
    #print full ANCOVA data to txt files
    np.savetxt(fileprefix + str(number_of_subjects) + '_' + str(len(TYPdata)) + 'TYP.txt', TYPdata, fmt = ["%.0f"]*2 + ["%.1f"] + ["%.6f"])
    np.savetxt(fileprefix + str(number_of_subjects) + '_' + str(len(ASDdata)) + 'ASD.txt', ASDdata, fmt = ["%.0f"]*2 + ["%.1f"] + ["%.6f"])
    #print IDs to txt files
    np.savetxt(subjectID_fileprefix + str(number_of_subjects) + '_' + str(len(TYPIDs)) + 'TYPids.txt', TYPIDs, fmt = ["%s"])
    np.savetxt(subjectID_fileprefix + str(number_of_subjects) + '_' + str(len(ASDIDs)) + 'ASDids.txt', ASDIDs, fmt = ["%s"])

def analyze_metrics(fileprefix,number_of_subjects,subjectID_fileprefix,outdir):
    #get data
    data, subIDs, metric = get_data(fileprefix,number_of_subjects,subjectID_fileprefix)
    #convert units mm^3 to cm^3
    data[:,metric['TCV']]=data[:,metric['TCV']]/1000
    #take max and min age and TCV
    maxage=np.max(data[:,metric['age']])
    minage=np.min(data[:,metric['age']])
    maxTCV=np.max(data[:,metric['TCV']])
    minTCV=np.min(data[:,metric['TCV']])

    #define ASD and TYP groups separately
    TYPyes=data[:,metric['diagnosis']]>0
    ASDyes=data[:,metric['diagnosis']]<0
    TYPdata=data[TYPyes]
    ASDdata=data[ASDyes]
    print "***********************************"
    print "******", sum(TYPyes), "TYPs,", sum(ASDyes), "ASDs ******"
    print "***********************************"

    #loop through each metric and print mean & stdev for each group, ttest of group comparison, plots
    for key in metric:
	asdkey=ASDdata[:,metric[key]]
	typkey=TYPdata[:,metric[key]]
	if key == 'gender':
    	    print ""
	    print "******", key, "******"
	    #print chi-square test for differences in gender distributions
	    print "%s percentage male: %6.3f" % ('ASD', np.mean(asdkey>0))
	    print "%s percentage male: %6.3f" % ('TYP', np.mean(typkey>0))
	    obs = np.array([[np.sum(asdkey>0), np.sum(asdkey<0)], [np.sum(typkey>0), np.sum(typkey<0)]])
	    chi2, p, dof, ex = stats.chi2_contingency(obs)
	    print 'chi2(%s) = %6.3f, p = %6.4f' %  (dof, chi2, p)
	elif key == 'diagnosis':
	    continue
	else:
	    print ""
	    print "******", key, "******"
	    #print mean and stdev of data
	    print "%s mean (sd): %4.2f (%4.2f)" % ('ASD', np.mean(asdkey), np.std(asdkey))
	    print "%s mean (sd): %4.2f (%4.2f)" % ('TYP', np.mean(typkey), np.std(typkey))
	    dof = sum(TYPyes)+sum(TYPyes==False)-2
	    t, p = stats.ttest_ind(asdkey, typkey)
	    print 'T(%s) = %6.3f, p = %6.4f' %  (dof, t, p)

    #define ASD-male, ASD-female, TYP-male, and TYP-female
    TYP_m_yes=TYPyes*(data[:,metric['gender']]>0)
    TYP_f_yes=TYPyes*(data[:,metric['gender']]<0)
    ASD_m_yes=ASDyes*(data[:,metric['gender']]>0)
    ASD_f_yes=ASDyes*(data[:,metric['gender']]<0)
    TYP_m_data=data[TYP_m_yes]
    TYP_f_data=data[TYP_f_yes]
    ASD_m_data=data[ASD_m_yes]
    ASD_f_data=data[ASD_f_yes]
    print ""
    print "***********************************"
    print "******", sum(TYP_m_yes), "TYP males,", sum(TYP_f_yes), "TYP females,", sum(ASD_m_yes), "ASD males,", sum(ASD_f_yes), "ASD females ******"
    print "***********************************"

    #name each group here
    groupdata = {
        'TYP males' : TYP_m_data,
        'TYP females' : TYP_f_data,
        'ASD males' : ASD_m_data,
        'ASD females' : ASD_f_data,
    }

    #loop through each metric and print mean & stdev for each group, ttest of group comparison
    for key in metric:
    	if key == 'diagnosis' or key == 'gender':
	    continue
	else:
	    print ""
	    print "******", key, "******"
	    #print mean and stdev of data for each group
	    for g in groupdata:
		meanval=np.mean(groupdata[g][:,metric[key]])
		stdval=np.std(groupdata[g][:,metric[key]])
		maxval=np.max(groupdata[g][:,metric[key]])
		submax=subIDs[data[:,metric[key]]==maxval]
		minval=np.min(groupdata[g][:,metric[key]])
		submin=subIDs[data[:,metric[key]]==minval]
		print "%s mean (sd): %4.2f (%4.2f), max is %s with %4.2f, min is %s with %4.2f" % (g, meanval, stdval, submax, maxval, submin, minval)
	    #compute t-tests of dignostic effect within each gender group
	    print "****** Tests of dignostic effect within each gender group ******"
	    print "T-test of TYP-male vs ASD-male"
	    dof = sum(TYP_m_yes)+sum(ASD_m_yes)-2
	    t, p = stats.ttest_ind(TYP_m_data[:,metric[key]], ASD_m_data[:,metric[key]])
	    print 'T(%s) = %6.3f, p = %6.4f' %  (dof, t, p)
	    print "T-test of TYP-female vs ASD-female"
	    dof = sum(TYP_f_yes)+sum(ASD_f_yes)-2
	    t, p = stats.ttest_ind(TYP_f_data[:,metric[key]], ASD_f_data[:,metric[key]])
	    print 'T(%s) = %6.3f, p = %6.4f' %  (dof, t, p)
	    #compute t-tests of gender effect within each diagnostic group
	    print "****** Tests of gender effect within each diagnostic group ******"
	    print "T-test of TYP-male vs TYP-female"
	    dof = sum(TYP_m_yes)+sum(TYP_f_yes)-2
	    t, p = stats.ttest_ind(TYP_m_data[:,metric[key]], TYP_f_data[:,metric[key]])
	    print 'T(%s) = %6.3f, p = %6.4f' %  (dof, t, p)
	    print "T-test of ASD-male vs ASD-female"
	    dof = sum(ASD_m_yes)+sum(ASD_f_yes)-2
	    t, p = stats.ttest_ind(ASD_m_data[:,metric[key]], ASD_f_data[:,metric[key]])
	    print 'T(%s) = %6.3f, p = %6.4f' %  (dof, t, p)

    #make awesome plots
    plt.figure()
    plt.subplot(211)
    plt.plot(TYP_m_data[:,metric['age']],TYP_m_data[:,metric['TCV']],'bo',TYP_f_data[:,metric['age']],TYP_f_data[:,metric['TCV']],'ro')
    plt.title('Typicals')
    plt.ylabel('TCV')
    plt.axis([minage-3,maxage+3,minTCV-50,maxTCV+50])
    plt.subplot(212)
    plt.plot(ASD_m_data[:,metric['age']],ASD_m_data[:,metric['TCV']],'bo',ASD_f_data[:,metric['age']],ASD_f_data[:,metric['TCV']],'ro')
    plt.title('ASD')
    plt.ylabel('TCV');plt.xlabel('age')
    plt.axis([minage-3,maxage+3,minTCV-50,maxTCV+50])
    plt.savefig(outdir + '/groupInfo_' + str(number_of_subjects) + '.png', dpi=300)
    plt.show()

#run the functions
mkdir_p(outdir)
analyze_metrics(fileprefix,number_of_subjects,subjectID_fileprefix,outdir)
make_gender_lists(fileprefix,filepreprefix,subjectID_fileprefix,number_of_subjects)
make_diag_lists(fileprefix,subjectID_fileprefix,number_of_subjects)
