# SRA ToolKit
Sequencing files are usually stored in databases such as NCBI and the ENA (European Nucleotide Archive). What if we want to download the sequenced and use them in our analysis? The SRA toolkit by NCBI is the solution in this case. SRA toolkit provides tools for downloading data, converting different formats of data into SRA dormat, and the vice versa, turning SRA data into different formats. 

The most frequently used tool from the toolkit would be fastq-dump that fetch SRA files and download them as fastq formats into our own directories. This page will provide information on how to install the SRA toolkit and use fastq-dump to download sequencing files for downstream analysis. 

# Installation
Please note the following instructions is for Ubuntu system only, and if you would like to install the toolkit for other systems, please download and install the correct package listed in this [page](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) for the system of choice to prevent configuration issues. 

1. Fetch the tar files from the NCBI: 
```
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
```

2. Extract the contents from the tar file:
```
tar -vxzf sratoolkit.tar.gz
```
- -v verbose shows the files tar works on while the command is running
- -x extract extract one or more items from an archive
- -z read or write compressed archives through bzip2 format
- -f specifies the file, follow by file names

3. we will append the path to the binaries to your PATH environment variable:
```
export PATH=$PATH:$PWD/sratoolkit.3.0.1-ubuntu64/bin
```

4. Before we configure the toolkit, we will verify the binaries will be found by the shell
```
which fastq-dump
```

5.  Now we can proceed to the configuration step followed by this command:
```
vdb-config i
```
These are the parameters we need to set up:
- Enable the "Remote Access" option on the Main screen.
- Go to cloud provider tab and accept to "report cloud instance identity".

6. Test that the toolkit if functional:
```
fastq-dump --stdout -X 2 SRR390728
```

the command should produce the following output
```
Read 2 spots for SRR390728
Written 2 spots for SRR390728
@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
```

# How to download multiple files at the same time
## The first thing we need is to create an accessing list, which is a file with all the sampleID you want to download. 
Here's one example that you can use: 

let nano this into a file called "accession_list.txt" 
```
SAMN07790141,SRR6178139,Camacho_Ortiz_2017
SAMN03002640,SRR1556545,PRJNA259188
```
The first column is the sampleID that we will feed into the fastq-dump to download, the second column is the runID that, and the third column is the name of the study just to keep an record. You can definetly add more information or less information for your purposes/own person habit but the sampleID is the required information. 

## Next we will feed the accession list to the fastq-dump to download multiple files at the same time
lets use nano to write this into a shell script called "fastq.sh"
```
while read line; do

# extract run number
SampleID=$(echo $line | cut -d',' -f1)
RunID=$(echo $line | cut -d',' -f2)
echo $SampleID $RunID

if [ ! -f "reads/${SampleID}/${RunID}_1.fastq.gz" ]; then
	mkdir -p reads/$SampleID # create a folder for each Biosample
	cd reads/$SampleID
	echo \-\-\-\> Downloading $SampleID
	fastq-dump --split-files --gzip $RunID # use fastq-dump to download its reads
	cd ../../
else
	echo \-\-\-\> $SampleID already downloaded
fi

done <scratch/accession_list.txt # move to next sample (ie line of file)
```

```
bash fastq.sh
```
