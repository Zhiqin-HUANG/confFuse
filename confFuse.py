#! /usr/bin/python
import os
import sys

################################################
#
# Name: confFuse
#
# Description: confFuse assigns a confidence score to  each putative fusion transcript from deFuse output
#			
# Author:  Zhiqin Huang, German Cancer Research Center (DKFZ)
# Date:	   2015-08-28
# Contact: z.huang@dkfz-heidelberg.de
#
################################################

###########
## required files: artefact list and gene annotation GTF
###########

# check input parameters
if len(sys.argv)!=5:
        print "Usage: python scripts.py defuse_folder output_folder artefact_list gencode.v19.annotation.gtf"
        sys.exit(0)

defuseFolder=sys.argv[1]
outputFolder=sys.argv[2]
artefactlistTable=sys.argv[3]
annotationFile=sys.argv[4]

# creat output folder if it does not exist
if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)
        print "create directory:" + outputFolder

# parameter weight, initial score = 10
# adjacent gene
adjacentW=0.5
# alternative splicing
altspliceW=4
# break point homology
breakpointsW=1
# exon boundary
exonboundariesW=1.5
# open reading frame
orfW=1
# read through, co-transcription
readThroughW=4
# interchromosomal translocation
interchrW=0.5
# split reads
splitW=0.5
# if fusion happens on downstream or utr3p of gene2
geneLocationW=4
# fusion in the same gene.
sameGeneW=4
# non-expression genes
zeroExpW=2

###################
# initiate score 10.
score=10
strict=10
# initiate list
artefactList=[]
occurrences=[]
# read artefact list table
for line in open(artefactlistTable,'r'):
	line=line.strip()
	artefactList.append(line)
# generating list of gene type based on annotation 
annotationDict=dict()
for line in open(annotationFile,'r'):
	lineParts=line.split("\"")
	if len(lineParts)>5:
		geneType=lineParts[5]
		geneID=lineParts[1].split(".")[0]
		annotationDict[geneID]=geneType

for filename in sorted(os.listdir(defuseFolder)):
	filePath=os.path.join(os.path.abspath(defuseFolder),filename)
	outputFile=os.path.join(os.path.abspath(outputFolder),filename + ".confFuse")
	output=open(outputFile,'w')

	with open(filePath,'r') as file_:
	        #calculate recurrence of fusion within each samples
		for line in file_:
                        line=line.strip()
                        lineParts=line.split("\t")
			gene1=lineParts[30]
                        gene2=lineParts[31]
                        if gene1 > gene2:
                                fusion=gene1 + "\t" + gene2
                        else:
                                fusion=gene2 + "\t" + gene1
			occurrences.append(fusion)		
		# sets the file's current position
		file_.seek(0)
		
		# calculate score
		firstLine=file_.readline()
		output.write(firstLine.strip() + "\t" + "geneType1" + "\t" + "geneType2" + "\t" + "occurrences" + "\t" + "inArtefactList" + "\t" + "adjustValue" + "\t" +  "confidence_score" + "\n")
		for line in file_:
			line=line.strip()
			lineParts=line.split("\t")
			split=lineParts[2]
			adj=lineParts[6]
			alts=lineParts[7]
			breakpoint=lineParts[11]
			exonb=lineParts[17]
			exp1=lineParts[18]
			exp2=lineParts[19]
			geneLocation1=lineParts[28]
                        geneLocation2=lineParts[29]
			interchr=lineParts[41]
			maxProportion=lineParts[47]
			numMultiMap=lineParts[50]
			orf=lineParts[52]
			readThough=lineParts[53]
			repeatProportion1=lineParts[54]
			repeatProportion1=lineParts[55]
			spanCount=lineParts[56]
			score=10	
			bonus=0
			strictBuff=0
			
			if adj=="Y":
				score=score-adjacentW
				strictBuff+=1
			if alts=="Y":
				score=score-altspliceW
			if int(breakpoint) >= 10:
				score=score-breakpointsW
			if exonb=="N":
				score=score-exonboundariesW
			if orf=="N":
				score=score-orfW
			if readThough=="Y":
				score=score-readThroughW
			if interchr=="N":
				score=score-interchrW	
			gene1=lineParts[30]
			gene2=lineParts[31]
			gene1ID=lineParts[20]
			gene2ID=lineParts[21]
			if gene1 > gene2:	
				fusion=gene1 + "\t" + gene2
			else:
				fusion=gene2 + "\t" + gene1
			if gene1==gene2:
				score=score-sameGeneW
			if int(occurrences.count(fusion)) >=2 and score<8 and float(maxProportion) < 0.5 and readThough=="N" and alts=="N":				
				if int(spanCount) - int(numMultiMap) + int(split) > 100 and float(numMultiMap)/float(spanCount)<0.9:
					score=score+1.5
					bonus+=1.5
					if geneLocation1=="coding":
						score=score+0.5
						bonus+=0.5
						if geneLocation2=="coding" or geneLocation2=="upstream":
							score=score+0.5
							bonus+=0.5
				elif int(spanCount) - int(numMultiMap) + int(split) > 40 and float(numMultiMap)/float(spanCount)<0.2:
					score=score+0.5
					bonus+=0.5
					if geneLocation1=="coding":
	                                	score=score+0.5
                                        	bonus+=0.5
                                                if geneLocation2=="coding" or geneLocation2=="upstream":
                                                        score=score+0.5
                                                        bonus+=0.5
			if geneLocation2=="downstream" or geneLocation2=="utr3p":
				score=score-geneLocationW
			if int(exp1)==0 or int(exp2)==0 and occurrences.count(fusion)<2 and int(spanCount) - int(numMultiMap) + int(split) <30:
				score=score-zeroExpW
                        if int(spanCount) - int(numMultiMap) < 5:
                                score=score-0.5
                                if int(numMultiMap) - int(spanCount) == 0:
                                        score=score-1
                                if int(split) < 5:
                                        score=score-splitW			
			if float(numMultiMap)/float(spanCount)>0.8 and int(numMultiMap)>5:
				score=score-0.5
			if float(maxProportion) > 0.9:
				score=score-1
			elif float(maxProportion) > 0.8:
				score=score-0.5
			gene1Type=annotationDict.get(gene1ID)
			gene2Type=annotationDict.get(gene2ID)
			if fusion in artefactList:
				score=score-6
				strictScore=score-strictBuff
				output.write(line + "\t" + str(gene1Type) + "\t" + str(gene2Type) + "\t" + str(occurrences.count(fusion)) + "\t" + "Y" + "\t" + str(bonus) + "\t" + str(strictScore) + "\n")
			else:
				strictScore=score-strictBuff
				output.write(line + "\t" + str(gene1Type) + "\t" + str(gene2Type) + "\t" + str(occurrences.count(fusion)) + "\t" + "N" + "\t" + str(bonus) + "\t" + str(strictScore) + "\n")
	occurrences=[]
	output.close()

