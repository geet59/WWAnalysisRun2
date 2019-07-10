import os
import glob
import os.path
searchfile = glob.glob('/uscms_data/d3/gchaudha/work_may/CMSSW_9_4_13/src/WWAnalysis/WWAnalysisRun2/OutPut_Logs/Logs_2019_07_04_06h00/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8_*stdout')
g= open('results.txt', 'w')
print(searchfile)
variable=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
variable_name = []
file_count = 0
for file_name in searchfile:
    print(file_name +"\n")
    f= open(file_name, 'r')
    count =0 
    for line in f:
     	if line.find("(",0,1) != -1 and line.find(")",1,5) != -1:
            g.write(line)
            print line,
            print "=====>  ",int(line.split(":")[1])
            if file_count==0: variable_name.append(line.split(":")[0])
            variable[count] += int(line.split(":")[1])
            print int(line.split(":")[1]), variable[count] 
            count+=1
            #for count in variable:
            #print(count) 
    file_count+=1
print "\n\n"
for length in range(0,len(variable_name)):
    print variable_name[length],variable[length]
f.close()

