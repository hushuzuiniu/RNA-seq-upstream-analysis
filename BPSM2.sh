#!/bin/bash
#Shu Hu

#fastqc
mkdir RNAseq
cp -r /localdisk/data/BPSM/AY21/fastq RNAseq/
cd RNAseq/fastq/
ls *gz |xargs fastqc -t 5
mkdir fastqc
mv  *fastqc.zip ./fastqc
mv *fastqc.html ./fastqc

#bulid index
cp /localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz .
gunzip TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz
hisat2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta TriTrypDB-46_TcongolenseIL3000_2019_Genome

gunzip *.gz

#maping and and converting "bam" format
for i in *_1.fq
do
i=${i/_1.fq/}
echo "hisat2 -p 20 -x TriTrypDB-46_TcongolenseIL3000_2019_Genome -1 ${i}_1.fq -2 ${i}_2.fq | samtools sort -@ 1 -o ${i}_sort.bam "
done > align_hisat2.sh
bash align_hisat2.sh

#bam index
for i in *_sort.bam
do
i=${i/_sort.bam/}
echo "samtools index  ${i}_sort.bam"
done > index.sh
bash index.sh

mkdir bam
mv *.bam bam
mv *.bam.bai bam
 
#bedtools
cd bam 
cp /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed TriTrypDB-46_TcongolenseIL3000_2019.bed
bedtools multicov -bams *.bam -bed TriTrypDB-46_TcongolenseIL3000_2019.bed > abc.txt

 cp /localdisk/data/BPSM/AY21/fastq/100k.fqfiles 100k.fqfiles

#calculate the means of each group 
function find_three_rows()
{

	let col=0
	let row=0
	let count=0
	let flag=0
	array=()

	for line in `cat 100k.fqfiles`
	do
	let col=${count}%7
	let row=${count}/7

	if [ ${col} == 0 ];
	then
		let temp_row=-1
		let flag=0
	fi

	if [[ ${col} == 1 ]]&&[[ ${line} == $1 ]];
	then
		temp_row=${row}
	fi

	if [[ ${row} == ${temp_row} ]]&&[[ ${col} == 3 ]];
	then
		if [ ${line} == $2 ];
		then
			let flag=${flag}+1
		fi
	fi

	if [[ ${row} == ${temp_row} ]]&&[[ ${col} == 4 ]];
	then
		if [ ${line} == $3 ];
		then
			let flag=${flag}+1
		fi
	fi

	if [ ${flag} == 2 ];
	then
		array[${#array[*]}]=${temp_row}
		let flag=0
	fi

	let count=${count}+1
	done

	echo ${array[*]}

}

function calculate_mean()
{
	line=$1

	let count=0
	for word in ${line}
	do
		data[${count}]=${word}
		let count=${count}+1
	done
	
	let offset=${count}-45
	
	let sum=0
	let index=0
	
	mean=()

	for element in $2
	do
		let mid_temp=${element}+${offset}-1
		let sum=${sum}+${data[${mid_temp}]}
		let index=${index}%3
		if [ ${index} == 2 ];
		then
			mean[${#mean[*]}]=$(echo ${sum} | awk '{printf ("%.2f",$1/3)}')
			let sum=0
		fi
		let index=${index}+1
	done
	
	echo ${mean[*]}
}

function process_raw_data()
{	
	row_array=$1
	
	cat "abc.txt" | while read line
	do
		mean=($(calculate_mean "${line}" "${row_array[*]}"))
		for element in ${mean[*]}
		do
			echo -n ${element} >> data.txt
			echo -n -e '\t' >> data.txt
		done
		echo >> data.txt
	done
}

#main function


Sample=("Clone1" "Clone2" "WT")
Time=("0" "24" "48")
Treatment=("Induced" "Uninduced")

let count=1

if [ -f "data.txt" ];
then
rm -r data.txt
fi
if [ ! -f "100k.fqfiles" ];
then
echo -e "\033[31mERROR:Can't find 100k.fqfiles\033[0m"
echo -e "\033[31mPlease put 100k.fqfiles and this script in the same path\033[0m"
exit
fi

if [ ! -f "abc.txt" ];
then
echo -e "\033[31mERROR:Can't find abc.txt\033[0m"
echo -e "\033[31mPlease put abc.txt and this script in the same path\033[0m"
exit
fi

echo "Running...please wait for several mins..."

for ((i=0;i< ${#Sample[*]};i ++))
do
	for ((j=0;j< ${#Time[*]};j ++))
	do
		for ((k=0;k< ${#Treatment[*]};k ++))
		do
		
		if [[ ${Time[j]} == "0" ]]&&[[ ${Treatment[k]} == "Induced" ]];
		then
			continue
		fi
		
		array=($(find_three_rows ${Sample[i]} ${Time[j]} ${Treatment[k]}))
		
		if [[ ${Sample[${i}]} == "Clone1" ]];
		then
			echo -n "C1_" >> data.txt
		elif [[ ${Sample[${i}]} == "Clone2" ]];
		then
			echo -n "C2_" >> data.txt
		elif [[ ${Sample[${i}]} == "WT" ]];
		then
			echo -n "WT_" >> data.txt
		fi
		
		echo -n ${Time[${j}]}"_" >> data.txt
		
		if [[ ${Treatment[${k}]} == "Induced" ]];
		then
			echo -n "I" >> data.txt
		elif [[ ${Treatment[${k}]} == "Uninduced" ]];
		then
			echo -n "U" >> data.txt
		fi
		
		echo -n -e '\t' >> data.txt
		
		let count=${count}+1
		ARRAY[${#ARRAY[*]}]=${array[0]}
		ARRAY[${#ARRAY[*]}]=${array[1]}
		ARRAY[${#ARRAY[*]}]=${array[2]}
		
		done
	done
done

echo >> data.txt

process_raw_data "${ARRAY[*]}"

#calculate the foldchanges and rank the results in descending order

if [ ! -f "data.txt" ];
then
echo -e "\033[31mERROR:Can't find data.txt\033[0m"
echo -e "\033[31mplease run process.sh first\033[0m"
exit
fi

echo "Running...please wait for several mins..."

if [ -f "output.txt" ];
then
rm -r output.txt
fi

if [ -f "temp.txt" ];
then
rm -r temp.txt
fi

if [ -f "temp_data.txt" ];
then
rm -r temp_data.txt
fi

if [ -f "output/temp.txt" ];
then
rm -r output/temp.txt
fi

if [ -d "output" ];
then
rm -r output
fi

awk 'NR>2{print line}{line=$0} END{print line}' data.txt > data_temp.txt

awk -F "\t" '{ print $1, $2, $3, $4, $5 }' abc.txt > temp.txt

mkdir output
mkdir output/compare_cells
mkdir output/compare_time
mkdir output/compare_treatment

awk '{printf ("%.2f\n",log(($1+1)/($11+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C1_0_U.txt

awk '{printf ("%.2f\n",log(($2+1)/($12+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C1_24_I.txt

awk '{printf ("%.2f\n",log(($3+1)/($13+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C1_24_U.txt

awk '{printf ("%.2f\n",log(($4+1)/($14+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C1_48_I.txt

awk '{printf ("%.2f\n",log(($5+1)/($15+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C1_48_U.txt

awk '{printf ("%.2f\n",log(($6+1)/($11+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C2_0_U.txt

awk '{printf ("%.2f\n",log(($7+1)/($12+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C2_24_I.txt

awk '{printf ("%.2f\n",log(($8+1)/($13+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C2_24_U.txt

awk '{printf ("%.2f\n",log(($9+1)/($14+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C2_48_I.txt

awk '{printf ("%.2f\n",log(($10+1)/($15+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_cells/C2_48_U.txt

#time

awk '{printf ("%.2f\n",log(($4+1)/($2+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_time/C1_I.txt

awk '{printf ("%.2f\n",log(($5+1)/($3+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_time/C1_U.txt

awk '{printf ("%.2f\n",log(($9+1)/($7+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_time/C2_I.txt

awk '{printf ("%.2f\n",log(($10+1)/($8+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_time/C2_U.txt

awk '{printf ("%.2f\n",log(($14+1)/($12+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_time/WT_I.txt

awk '{printf ("%.2f\n",log(($15+1)/($13+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_time/WT_U.txt

#treatment
awk '{printf ("%.2f\n",log(($2+1)/($3+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_treatment/C1_24.txt

awk '{printf ("%.2f\n",log(($4+1)/($5+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_treatment/C1_48.txt

awk '{printf ("%.2f\n",log(($7+1)/($8+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_treatment/C2_24.txt

awk '{printf ("%.2f\n",log(($9+1)/($10+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_treatment/C2_48.txt

awk '{printf ("%.2f\n",log(($12+1)/($13+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_treatment/WT_24.txt

awk '{printf ("%.2f\n",log(($14+1)/($15+1))/log(10))}' data_temp.txt > temp_data.txt
paste temp.txt temp_data.txt > output/temp.txt
awk '{print $NF,$0}' output/temp.txt | sort -n -r| cut -f2- -d' ' > output/compare_treatment/WT_48.txt


rm temp_data.txt
rm temp.txt
rm output/temp.txt
rm data_temp.txt

echo "Complete! The result is output folder"