#!/bin/bash
set -e;

function usage {
    echo -e "usage : parseGTFAnnotation.sh -i GTF [-g GENOME_SIZE] [-o OUTPUT_PATH] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "parseGencodeAnnotation.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i GTF: input GTF file from Gencode"
    echo "   [-g GENOME_SIZE]: chromosome sizes file"
    echo "   [-o OUTPUT_PATH]: output path"
    echo "   [-h|--help]: help"
    exit;
}

while getopts "i:g:o:h" OPT
do
    case $OPT in
        i) INPUT=$OPTARG;;
	g) GENOME_SIZE=$OPTARG;;
	o) ODIR=$OPTARG;;
	h) help ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    usage
	     exit 1
	     ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    usage
	    exit 1
	    ;;
    esac
done


if [[ -z $INPUT ]]; then
    usage
    exit
fi

if [[ -z $ODIR ]]; then
    ODIR="./"
fi

if [[ ${INPUT} =~ ".gff" ]]; then
    echo -e "GFF format detected ..."
    FORMAT="gff"
elif [[ ${INPUT} =~ ".gtf" ]]; then
    echo -e "GTF format detected ..."
    FORMAT="gtf"
else
    exit 1
fi

## Define outputs
mkdir -p ${ODIR}
if [[ $FORMAT == "gtf" ]]; then
    OUTPUT_EXON=${ODIR}/$(basename ${INPUT} | sed -e 's/.gtf$/_exon.bed/')
    OUTPUT_GENE=${ODIR}/$(basename ${INPUT} | sed -e 's/.gtf$/_gene.bed/')
    OUTPUT_PROMOTER=${ODIR}/$(basename ${INPUT} | sed -e 's/.gtf$/_promoter_2000.bed/')
    OUTPUT_INTRON=${ODIR}/$(basename ${INPUT} | sed -e 's/.gtf$/_intron.bed/')
    OUTPUT_INTER=${ODIR}/$(basename ${INPUT} | sed -e 's/.gtf$/_inter.bed/')
    OUTPUT_GENEPRED=${ODIR}/$(basename ${INPUT} | sed -e 's/.gtf$/.genePred/')
    OUTPUT_BED12=${ODIR}/$(basename ${INPUT} | sed -e 's/.gtf$/.bed12/')
else
    OUTPUT_EXON=${ODIR}/$(basename ${INPUT} | sed -e 's/.gff[0-9]$/_exon.bed/')
    OUTPUT_GENE=${ODIR}/$(basename ${INPUT} | sed -e 's/.gff[0-9]$/_gene.bed/')
    OUTPUT_PROMOTER=${ODIR}/$(basename ${INPUT} | sed -e 's/.gff[0-9]$/_promoter_2000.bed/')
    OUTPUT_INTRON=${ODIR}/$(basename ${INPUT} | sed -e 's/.gff[0-9]$/_intron.bed/')
    OUTPUT_INTER=${ODIR}/$(basename ${INPUT} | sed -e 's/.gff[0-9]$/_inter.bed/')
    OUTPUT_GENEPRED=${ODIR}/$(basename ${INPUT} | sed -e 's/.gff[0-9]$/.genePred/')
    OUTPUT_BED12=${ODIR}/$(basename ${INPUT} | sed -e 's/.gff[0-9]$/.bed12/')
fi


## Get genes
echo -e "Extract genic regions ..."
if [[ ${FORMAT} == "gff" ]]; then
    extractGeneFromGTF.R ${INPUT} ${FORMAT} ${OUTPUT_GENE}_tmp 2> extractGeneFromGTF.Rout
    sort -k1,1V -k2,2n -k3,3n ${OUTPUT_GENE}_tmp > ${OUTPUT_GENE}
    rm -f ${OUTPUT_GENE}_tmp
elif [[ ${FORMAT} == "gtf" ]]; then
    if [[ $(head -500 ${INPUT} | awk '$3=="gene"{print}' | wc -l) -eq "0" ]]; then
	echo -e "Warning: No gene information detected ... swith to R extraction ..."
        extractGeneFromGTF.R ${INPUT} ${FORMAT} ${OUTPUT_GENE}_tmp 2> extractGeneFromGTF.Rout
	sort -k1,1V -k2,2n -k3,3n ${OUTPUT_GENE}_tmp > ${OUTPUT_GENE}
	rm -f ${OUTPUT_GENE}_tmp
    else
	## Deal with special case for FlyBase/SGD
	awk -F"\t" 'BEGIN{OFS="\t"} $3=="gene"{split($9,annot,";"); for (i=1;i<=length(annot);i++){split(annot[i],b," "); hannot[b[1]]=b[2]; gsub("\"","",hannot[b[1]])}; if ($2=="FlyBase"){n=hannot["gene_symbol"]} else if ($2=="SGD"){n=hannot["gene_id"]"|"hannot["gene_name"]; gsub("\\|$","",n)} else{n=hannot["gene_name"]}; print $1,$4-1,$5,n,"0",$7; delete hannot}' ${INPUT} | sort -k1,1V -k2,2n -k3,3n > ${OUTPUT_GENE}
    fi
fi

## Get promoters
echo -e "Extract promoter regions ..."
awk -F"\t" -v win=2000 'BEGIN{OFS="\t"} $6=="+"{s=$2-win;e=$2+win;if(s<0){s=0}; print $1,s,e,$4,$5,$6} $6=="-"{print $1,$3-win,$3+win,$4,$5,$6}' ${OUTPUT_GENE} > ${OUTPUT_PROMOTER}

## Get exonic regions
echo -e "Extract exonic regions ..."
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' ${INPUT} | sort -k1,1V -k2,2n -k3,3n > ${OUTPUT_EXON}
mergeBed -i ${OUTPUT_EXON} > ${OUTPUT_EXON}_temp
mv ${OUTPUT_EXON}_temp ${OUTPUT_EXON}

## Get Intronic regions
echo -e "Extract intronic regions ..."
subtractBed -a ${OUTPUT_GENE} -b ${OUTPUT_EXON} > ${OUTPUT_INTRON}

if [[ ! -z ${GENOME_SIZE} ]]; then
    ## Get Intergenic regions
    ## BEDtools 2.25.0: Bug issue with end > start
    ## BEDtools 2.27.1: Bug issue with end = start
    ## BEDtools 2.27.1: crash with all chromosomes on mm10
    echo -e "Extract intergenic regions ..."
    sort -k1,1V ${GENOME_SIZE} > genome_size.temp
    complementBed -i ${OUTPUT_GENE} -g genome_size.temp | awk '$2<$3 {print}' > ${OUTPUT_INTER}
    rm genome_size.temp
fi

## Generate bed12 file from gtf
## Based on UCSC tools
echo -e "Generate bed12 file ..."
if [[ $FORMAT == "gtf" ]]; then
    gtfToGenePred -allErrors -ignoreGroupsWithoutExons ${INPUT} ${OUTPUT_GENEPRED} 2> ${ODIR}/genepred.log
    genePredToBed ${OUTPUT_GENEPRED} ${OUTPUT_BED12}
    rm ${OUTPUT_GENEPRED}
elif [[ $FORMAT == "gff" ]]; then
    gff3ToGenePred ${INPUT} ${OUTPUT_GENEPRED} 2> ${ODIR}/genepred.log
    genePredToBed ${OUTPUT_GENEPRED} ${OUTPUT_BED12}
    rm ${OUTPUT_GENEPRED}
fi
