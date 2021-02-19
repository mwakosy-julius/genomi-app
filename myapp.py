import streamlit as st
import pandas as pd
import altair as alt
from PIL import Image

st.title('WELCOME TO MY PETITE GENOMI APP')

#image = Image.open('Untitled.jpg')
#st.image(image, use_column_width=True)

st.header('Enter Sequence In CAPS')

sequence = st.text_area("", height=25)

DNA = st.button('DNA')
RNA = st.button('RNA')

def dna_dna():
    replica = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'T'
        elif nucleotide == 'T':
            nucleotide = 'A'
        elif nucleotide == 'C':
            nucleotide = 'G'
        elif nucleotide == 'G':
            nucleotide = 'C'
        else:
            break
        replica+=nucleotide
    return replica

def transcription():
    transcript = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'U'
        elif nucleotide == 'T':
            nucleotide = 'A'
        elif nucleotide == 'C':
            nucleotide = 'G'
        elif nucleotide == 'G':
            nucleotide = 'C'
        else:
            break
        transcript+=nucleotide
    return transcript

def central():
    sequence2 = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'U'
        elif nucleotide == 'T':
            nucleotide = 'A'
        elif nucleotide == 'C':
            nucleotide = 'G'
        elif nucleotide == 'G':
            nucleotide = 'C'
        sequence2 += nucleotide

    RNA_sequence = []
    n = 0
    k = 1
    for seq in sequence2:
        RNA_sequence.append(sequence2[n]+sequence2[n+1]+sequence2[n+2])
        if len(sequence2)//3 > k:
            n += 3
            k += 1
        else:
            break

    amino_sequence = ''
    for codon in RNA_sequence:
        if codon == 'UUU' or codon == 'UUC':
            codon = 'Phe-'
        elif codon == 'UUA' or codon == 'UUG' or codon == 'CUU' or codon == 'CUC' or codon == 'CUA' or codon == 'CUG':
            codon = 'Leu-'
        elif codon == 'AUU' or codon == 'AUC' or codon == 'AUA':
            codon = 'Ile-'
        elif codon == 'AUG':
            codon = 'Met-'
        elif codon == 'GUU' or codon == 'GUC' or codon == 'GUA' or codon == 'GUG':
            codon = 'Val-'
        elif codon == 'UCU' or codon == 'UCC' or codon == 'UCA' or codon == 'UCG' or codon == 'AGU' or codon == 'AGC':
            codon = 'Ser-'
        elif codon == 'CCU' or codon == 'CCC' or codon == 'CCA' or codon == 'CCG':
            codon = 'Pro-'
        elif codon == 'ACU' or codon == 'ACC' or codon == 'ACA' or codon == 'ACG':
            codon = 'Thr-'
        elif codon == 'GCU' or codon == 'GCC' or codon == 'GCA' or codon == 'GCG':
            codon = 'Ala-'
        elif codon == 'UAU' or codon == 'UAC':
            codon = 'Tyr-'
        elif codon == 'UAA' or codon == 'UAG' or codon == 'UGA':
            codon = 'STOP'
            break
        elif codon == 'UAU' or codon == 'UAC':
            codon = 'Tyr-'
        elif codon == 'CAU' or codon == 'CAC':
            codon = 'His-'
        elif codon == 'CAA' or codon == 'CAG':
            codon = 'Gln-'
        elif codon == 'AAU' or codon == 'AAC':
            codon = 'Asn-'
        elif codon == 'AAA' or codon == 'AAG':
            codon = 'Lys-'
        elif codon == 'GAU' or codon == 'GAC':
            codon = 'Asp-'
        elif codon == 'GAA' or codon == 'GAG':
            codon = 'Glu-'
        elif codon == 'UGU' or codon == 'UGC':
            codon = 'Cys-'
        elif codon == 'UGG':
            codon = 'Trp-'
        elif codon == 'CGU' or codon == 'CGC' or codon == 'CGA' or codon == 'CGG' or codon == 'AGA' or codon == 'AGG':
            codon = 'Arg-'
        elif codon == 'GGU' or codon == 'GGC' or codon == 'GGA' or codon == 'GGG':
            codon = 'Gly-'
        amino_sequence+=codon
    return amino_sequence

def percentaged():
    d = dict([
    ('A ', sequence.count('A')),
    ('G ', sequence.count('G')),
    ('C ', sequence.count('C')),
    ('T ', sequence.count('T'))
    ])
    return d

def percentager():
    r = dict([
    ('A ', sequence.count('A')),
    ('G ', sequence.count('G')),
    ('C ', sequence.count('C')),
    ('U ', sequence.count('U'))
    ])
    return r

def translation():
    RNA_sequence = []
    n = 0
    k = 1
    for seq in sequence:
        RNA_sequence.append(sequence[n]+sequence[n+1]+sequence[n+2])
        if len(sequence)//3 > k:
            n += 3
            k += 1
        else:
            break

    amino_sequence = ''
    for codon in RNA_sequence:
        if codon == 'UUU' or codon == 'UUC':
            codon = 'Phe-'
        elif codon == 'UUA' or codon == 'UUG' or codon == 'CUU' or codon == 'CUC' or codon == 'CUA' or codon == 'CUG':
            codon = 'Leu-'
        elif codon == 'AUU' or codon == 'AUC' or codon == 'AUA':
            codon = 'Ile-'
        elif codon == 'AUG':
            codon = 'Met-'
        elif codon == 'GUU' or codon == 'GUC' or codon == 'GUA' or codon == 'GUG':
            codon = 'Val-'
        elif codon == 'UCU' or codon == 'UCC' or codon == 'UCA' or codon == 'UCG' or codon == 'AGU' or codon == 'AGC':
            codon = 'Ser-'
        elif codon == 'CCU' or codon == 'CCC' or codon == 'CCA' or codon == 'CCG':
            codon = 'Pro-'
        elif codon == 'ACU' or codon == 'ACC' or codon == 'ACA' or codon == 'ACG':
            codon = 'Thr-'
        elif codon == 'GCU' or codon == 'GCC' or codon == 'GCA' or codon == 'GCG':
            codon = 'Ala-'
        elif codon == 'UAU' or codon == 'UAC':
            codon = 'Tyr-'
        elif codon == 'UAA' or codon == 'UAG' or codon == 'UGA':
            codon = 'STOP'
            break
        elif codon == 'UAU' or codon == 'UAC':
            codon = 'Tyr-'
        elif codon == 'CAU' or codon == 'CAC':
            codon = 'His-'
        elif codon == 'CAA' or codon == 'CAG':
            codon = 'Gln-'
        elif codon == 'AAU' or codon == 'AAC':
            codon = 'Asn-'
        elif codon == 'AAA' or codon == 'AAG':
            codon = 'Lys-'
        elif codon == 'GAU' or codon == 'GAC':
            codon = 'Asp-'
        elif codon == 'GAA' or codon == 'GAG':
            codon = 'Glu-'
        elif codon == 'UGU' or codon == 'UGC':
            codon = 'Cys-'
        elif codon == 'UGG':
            codon = 'Trp-'
        elif codon == 'CGU' or codon == 'CGC' or codon == 'CGA' or codon == 'CGG' or codon == 'AGA' or codon == 'AGG':
            codon = 'Arg-'
        elif codon == 'GGU' or codon == 'GGC' or codon == 'GGA' or codon == 'GGG':
            codon = 'Gly-'
        amino_sequence+=codon
    return amino_sequence

def rna_rna():
    replica = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'U'
        elif nucleotide == 'U':
            nucleotide = 'A'
        elif nucleotide == 'C':
            nucleotide = 'G'
        elif nucleotide == 'G':
            nucleotide = 'C'
        replica += nucleotide
    return replica

def rna_dna():
    reverse = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'T'
        elif nucleotide == 'U':
            nucleotide = 'A'
        elif nucleotide == 'C':
            nucleotide = 'G'
        elif nucleotide == 'G':
            nucleotide = 'C'
        reverse += nucleotide
    return reverse



replication = dna_dna()
transcription = transcription()
aminoAcids = central()
percentages = percentaged()

rnaReplication = rna_rna()
reverse_transcription = rna_dna()
translation = translation()
percentager = percentager()

def valid():
    chance = 0
    for nucleotide in sequence:
        if nucleotide == 'A' or nucleotide == 'T' or nucleotide == 'G' or nucleotide == 'C':
            chance+=1
        else:
            chance-=1
    if chance == len(sequence) and chance > 0:
        st.success('Your DNA is dynamite')
        st.write('DNA sequence: ',  sequence)
        st.write('DNA replicon: ', replication)
        st.write('DNA transript: ', transcription)
        st.write('Amino acids: ', aminoAcids)
        st.write('Nucleotide Count: ', percentages)

    if chance < len(sequence):
        st.error('Ooops, What kind of DNA is this!')
    if chance == 0:
        st.warning("Hey don't just look there, Enter something")
        
def valir():
    chance = 0
    for nucleotide in sequence:
        if nucleotide == 'A' or nucleotide == 'U' or nucleotide == 'G' or nucleotide == 'C':
            chance+=1
        else:
            chance-=1
    if chance == len(sequence) and chance > 0:
        st.success('Your RNA is dynamite')
        st.write('RNA sequence: ',  sequence)
        st.write('RNA replicon: ', rnaReplication)
        st.write('DNA transript: ', reverse_transcription)
        st.write('Amino acids: ', translation)
        st.write('Nucleotide Count: ', percentager)
        
    if chance < len(sequence):
        st.error('Ooops, What kind of RNA is this!')
    if chance == 0:
        st.warning("Hey don't just look there, Enter something")

if DNA:
    valid() 
if RNA:
    valir()

