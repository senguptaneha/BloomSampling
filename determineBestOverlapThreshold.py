import sys
import os

MValues = [100000,1000000];
MStrings = ['100K','1M'];

nValues = [100,1000,10000,50000];
nStrings = ['100','1K','10K','50K'];
#mValues = [[62000,64000,80000,88000,100000,128000],[620000,640000,760000,800000,1000000,1280000]];
mValues = [[6313,7248,8422,10110,13294,25000],
	   [28465,32808,38259,46000,60870,100000],
	   [123175,143444,168630,204938,273405,300000],
	   [319413,378200,451743,557011,755093,900000]];

levelValues = [[9,9,10,12,8,10],
	       [10,10,9,9,9,9],
	       [11,8,8,8,7,7],
	       [13,13,13,6,6,6]];

for i in range(1,2):
	M = MValues[i];
	for j in range(3,4):
		for z in range(6):
			m = mValues[j][z];
			lT = levelValues[j][z];
			minTime = sys.float_info.max;
			bestOverlap = 0;
			bestIntersection = 0;
			bestMembership = 0;
			k = 0;
			while (k < 2):
				k = k + 0.01;
				f = os.popen('./reconstructBT '+`M`+' '+`nValues[j]`+' ../SimulatedDatasets/'+MStrings[i]+'/Skewed'+nStrings[j]+'_1 3 '+`m`+' '+`k`+' '+`lT`);
				l = f.read();
				f.close();
				print(`k`+','+`nValues[j]`+':'+l);

				#print(`k`+':'+l);
				precisionTok = (l.split('PRECISION ')[1]).split(' ')[0];
				if (precisionTok.find('inf')!=-1):
					break;
				precision = float(precisionTok);
				if (precision >1): continue;
				membership = int((l.split('MEMBERSHIP=')[1]).split(' ')[0]);
				if (membership <= 0):
					continue;
				time = float((l.split('TIME ')[1]).split(' ')[0]);
				if (time<minTime):
					minTime = time;
					bestOverlap = k;
					bestIntersection = int((l.split('INTERSECTIONS=')[1]).split(',')[0]);
					bestMembership = membership;
				
			print(`M`+'\t'+`m`+'\t'+`nValues[j]`+'\t'+`bestOverlap`+'\t'+`bestIntersection`+'\t'+`bestMembership`+'\t'+`minTime`);
					
							
