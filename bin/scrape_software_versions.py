#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'Pipeline': ['v_main.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'Bowtie2': ['v_bowtie2.txt', r"version (\S+)"],
    'STAR': ['v_star.txt', r"(\S+)"],
    'HISAT2': ['v_hisat2.txt', r"version (\S+)"],
    'TopHat2': ['v_tophat2.txt', r"TopHat v(\S+)"],
    'Picard MarkDuplicates': ['v_markduplicates.txt', r"([\d\.]+)-SNAPSHOT"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'featureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'deeptools': ['v_deeptools.txt', r"plotFingerprint (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['Pipeline'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie2'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = False
results['HISAT2'] = False
results['TopHat2'] = '<span style="color:#999999;\">N/A</span>'
results['Picard MarkDuplicates'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['featureCounts'] = '<span style="color:#999999;\">N/A</span>'
results['deeptools'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = match.group(1)

# Strip STAR or HiSAT2
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'Software Versions'
section_href: 'https://gitlab.curie.fr/rnaseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")
