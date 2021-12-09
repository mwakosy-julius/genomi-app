import streamlit as st
import pandas as pd
import altair as alt
from PIL import Image

st.title('WELCOME TO MY GENOMI APP')

#image = Image.open('chroma.jpg')
#st.image(image, use_column_width=True)

col1, col2 = st.beta_columns(2)

st.header('Enter Sequence')

sequence = st.text_area("", height=25)
sequence = sequence.upper()
nott = '\n'
sequence = sequence.replace(nott, '')
sequencer = ' '.join(sequence)

DNA = st.button('DNA')
RNA = st.button('RNA')

def dna_dna():
    replica = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'T '
        elif nucleotide == 'T':
            nucleotide = 'A '
        elif nucleotide == 'C':
            nucleotide = 'G '
        elif nucleotide == 'G':
            nucleotide = 'C '
        else:
            break
        replica+=nucleotide
    
    return replica

def transcription():
    transcript = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'U '
        elif nucleotide == 'T':
            nucleotide = 'A '
        elif nucleotide == 'C':
            nucleotide = 'G '
        elif nucleotide == 'G':
            nucleotide = 'C '
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

aminoAcids = central()

def gc():
    gc = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    content = ''
    try:
        content = str(int(gc/total*100))+'%'
    except ZeroDivisionError:
        pass
    return content

def total_count():
    duuh = str(len(sequence))
    return duuh

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

def aminod(seq):
    acid = dict([
    ('Phe ', seq.count('Phe')),
    ('Leu ', seq.count('Leu')),
    ('Ile ', seq.count('Ile')),
    ('Met ', seq.count('Met')),
    ('Val ', seq.count('Val')),
    ('Ser ', seq.count('Ser')),
    ('Pro ', seq.count('Pro')),
    ('Thr ', seq.count('Thr')),
    ('Ala ', seq.count('Ala')),
    ('Tyr ', seq.count('Tyr')),
    ('His ', seq.count('HIs')),
    ('Gln ', seq.count('Gln')),
    ('Asn ', seq.count('Asn')),
    ('Lys ', seq.count('Lys')),
    ('Asp ', seq.count('Asp')),
    ('Glu ', seq.count('Glu')),
    ('Cys ', seq.count('Cys')),
    ('Trp ', seq.count('Trp')),
    ('Arg ', seq.count('Arg')),
    ('Gly ', seq.count('Gly'))
    ])
    return acid
number = aminod(aminoAcids)

def graphd():
    X = percentaged()
    X_label = list(X)
    X_values = list(X.values())

    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': 'nucleotide'})
    return df

def graphr():
    X = percentager()
    X_label = list(X)
    X_values = list(X.values())

    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': 'nucleotide'})
    return df

def chartd():
    X = percentaged()
    X_label = list(X)
    X_values = list(X.values())

    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': 'nucleotide'})
    

    p = alt.Chart(df).mark_bar().encode(
        x = 'nucleotide',
        y = 'count'
    )

    p = p.properties(
        width = alt.Step(80)
    )
    return p

def chartr():
    X = percentager()
    X_label = list(X)
    X_values = list(X.values())

    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': 'nucleotide'})
    

    p = alt.Chart(df).mark_bar().encode(
        x = 'nucleotide',
        y = 'count'
    )

    p = p.properties(
        width = alt.Step(80)
    )
    return p

def charta(num):
    X = num
    X_label = list(X)
    X_values = list(X.values())

    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': 'amino-acid'})
    

    p = alt.Chart(df).mark_bar().encode(
        x = 'amino-acid',
        y = 'count'
    )

    p = p.properties(
        width = alt.Step(20)
    )
    return p
chacha = charta(number)

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

aminoAcid = translation()
numbers = aminod(aminoAcid)
chachas = charta(numbers)

def rna_rna():
    replica = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'U '
        elif nucleotide == 'U':
            nucleotide = 'A '
        elif nucleotide == 'C':
            nucleotide = 'G '
        elif nucleotide == 'G':
            nucleotide = 'C '
        replica += nucleotide
    return replica

def rna_dna():
    reverse = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            nucleotide = 'T '
        elif nucleotide == 'U':
            nucleotide = 'A '
        elif nucleotide == 'C':
            nucleotide = 'G '
        elif nucleotide == 'G':
            nucleotide = 'C '
        reverse += nucleotide
    return reverse

replication = dna_dna()
transcription = transcription()
gc_content = gc()
total_count = total_count()
percentages1 = graphd()
percentages2 = chartd()
amino_chart = charta(number)

rnaReplication = rna_rna()
reverse_transcription = rna_dna()
translation = translation()
percentagesa = graphr()
percentagesb = chartr()

def valid():
    chance = 0
    for nucleotide in sequence:
        if nucleotide == 'A' or nucleotide == 'T' or nucleotide == 'G' or nucleotide == 'C':
            chance+=1
        elif nucleotide == ' ':
            nucleotide = ''
            chance+=1
        else:
            st.error('Invalid input in nucleotide number '+ str(sequence.index(nucleotide)+1))
            break
            chance-=1
    
    if len(sequence) == 0:
        st.warning("Hey don't just look there, Enter something")

    if chance == len(sequence) and chance > 0:
        st.success('Your DNA is fantastic')
        st.write('dna sequence   \n',  sequencer)
        st.write('dna replicon   \n', replication)
        st.write('dna transcript   \n ', transcription)
        st.write('AMINO ACIDS   \n', aminoAcids)
        st.write(percentages1)
        st.write(percentages2)
        st.write('The G-C content in the sequence is ', gc_content)
        st.write('The total number of nucleotides is ', total_count)
        st.write(chacha)

    
        
def valir():
    chance = 0
    for nucleotide in sequence:
        if nucleotide == 'A' or nucleotide == 'U' or nucleotide == 'G' or nucleotide == 'C':
            chance+=1
        elif nucleotide == ' ':
            nucleotide = ''
            chance+=1
        else:
            st.error('Invalid input in nucleotide number '+ str(sequence.index(nucleotide)+1))
            break
            chance-=1
    
    if len(sequence) == 0:
        st.warning("Hey don't just look there, Enter something")

    if chance == len(sequence) and chance > 0:
        st.success('Your RNA is fantastic')
        st.write('rna sequence   \n',  sequencer)
        st.write('rna replicon   \n', rnaReplication)
        st.write('reverse transript   \n', reverse_transcription)
        st.write('Amino acids   \n', translation)
        st.write(percentagesa)
        st.write(percentagesb)
        st.write('The G-C content in the sequence is ', gc_content)
        st.write('The total number of nucleotides is ', total_count)
        st.write(chachas)
        
if DNA:
    valid() 
if RNA:
    valir()







