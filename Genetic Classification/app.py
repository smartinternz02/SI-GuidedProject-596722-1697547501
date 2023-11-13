from flask import Flask, render_template, request
import pickle
import numpy as np
import pandas as pd

app = Flask(__name__)

model = pickle.load(open('model.pkl', 'rb'))
scaler = pickle.load(open('Scaler.pkl', 'rb'))
encoder_1 = pickle.load(open('Encoder_1.pkl', 'rb'))
encoder_2 = pickle.load(open('Encoder_2.pkl', 'rb'))
encoder_3 = pickle.load(open('Encoder_3.pkl', 'rb'))
encoder_4 = pickle.load(open('Encoder_4.pkl', 'rb'))
encoder_5 = pickle.load(open('Encoder_5.pkl', 'rb'))
encoder_6 = pickle.load(open('Encoder_6.pkl', 'rb'))
encoder_7 = pickle.load(open('Encoder_7.pkl', 'rb'))
encoder_8 = pickle.load(open('Encoder_8.pkl', 'rb'))
encoder_9 = pickle.load(open('Encoder_9.pkl', 'rb'))
encoder_10 = pickle.load(open('Encoder_10.pkl', 'rb'))

@app.route('/')
def start():
    return render_template('index.html')


@app.route('/login', methods=['POST'])
def login():
    chrom = request.form["ch"] or np.NaN
    pos = request.form["ps"] or np.NaN
    ref = request.form["rf"] or np.NaN
    alt = request.form["at"] or np.NaN
    af_esp = request.form["ap"] or np.NaN
    af_exac = request.form["ac"] or np.NaN
    af_tgp = request.form["ag"] or np.NaN
    clndisdb = request.form["cb"] or np.NaN
    clndisdbincl = request.form["cl"] or np.NaN
    clndn = request.form["cn"] or np.NaN
    clndnincl = request.form["cc"] or np.NaN
    clnhgvs = request.form["cs"] or np.NaN
    clnsigincl = request.form["ci"] or np.NaN
    clnvc = request.form["cv"] or np.NaN
    clnvi = request.form["li"] or np.NaN
    mc = request.form["mc"] or np.NaN
    origin = request.form["on"] or np.NaN
    ssr = request.form["sr"] or np.NaN
    allele = request.form["ae"] or np.NaN
    consequence = request.form["cq"] or np.NaN
    impact = request.form["im"] or np.NaN
    symbol = request.form["sl"] or np.NaN
    feature_type = request.form["fet_type"] or np.NaN
    feature = request.form["fet"] or np.NaN
    biotype = request.form["bt"] or np.NaN
    exon = request.form["ex"] or np.NaN
    intron = request.form["int"] or np.NaN
    cdna_position = request.form["cDNA"] or np.NaN
    cds_position = request.form["cds"] or np.NaN
    protein_position = request.form["pp"] or np.NaN
    amino_acids = request.form["aa"] or np.NaN
    codons = request.form["cd"] or np.NaN
    distance = request.form["ds"] or np.NaN
    strand = request.form["sd"] or np.NaN
    bam_edit = request.form["be"] or np.NaN
    sift = request.form["st"] or np.NaN
    polyphen = request.form["php"] or np.NaN
    motif_name = request.form["mn"] or np.NaN
    motif_pos = request.form["mp"] or np.NaN
    high_inf_pos = request.form["hip"] or np.NaN
    motif_score_change = request.form["msc"] or np.NaN
    loftool = request.form["lft"] or np.NaN
    cadd_phred = request.form["cp"] or np.NaN
    cadd_raw = request.form["cr"] or np.NaN
    blosum62 = request.form["bls"] or np.NaN

    if ref not in ['A', 'C', 'G', 'T']:
        ref = 'Other'
    if alt not in ['A', 'C', 'G', 'T']:
        alt = 'Other'
    if allele not in ['A', 'C', 'G', 'T']:
        allele = 'Other'

    if clnvc in ['Insertion', 'Inversion', 'Microsatellite']:
        clnvc = 'Other'

    if pd.isnull(clnvi):
        clnvi = 0
    else:            
        clnvi = 1

    if pd.isnull(intron):
        intron = 0
    else:            
        intron = 1

    if pd.isnull(bam_edit):
        bam_edit = 0
    else:            
        bam_edit = 1

    if pd.isnull(sift):
        sift = 0
    else:            
        sift = 1

    if pd.isnull(polyphen):
        polyphen = 0
    else:            
        polyphen = 1

    if pd.isnull(blosum62):
        blosum62 = 0
    else:            
        blosum62 = 1    

    if origin != '1':
        origin = 'Other'
    
    if symbol not in ['TTN', 'BRCA2', 'ATM', 'APC', 'BRCA1', 'MSH6', 'LDLR', 'PALB2', 'NF1', 'TSC2', 'BRIP1', 'PMS2', 'MSH2', 'POLE', 'CDH1', 'CHEK2', 'BARD1', 'SMARCA4', 'MYBPC3', 'RAD50', 'SYNE1', 'TP53', 'POLD1', 'NBN', 'MLH1', 'MUTYH', 'PLEC', 'STK11', 'MYH7', 'NEB', 'DMD', 'RYR2', 'RYR1', 'FBN1', 'DICER1', 'SCN5A', 'COL6A3', 'PTEN', 'RAD51C', 'RAD51D', 'DSP', 'RET', 'USH2A', 'TSC1', 'DYSF', 'CDH23', 'PTCH1', 'COL5A1',
                      'SYNE2', 'BMPR1A', 'DNAH11', 'FBN2', 'COL6A2', 'GPR98', 'MYH11', 'MYH6', 'COL6A1', 'CFTR', 'GAA', 'POLG', 'PKHD1', 'DNAH5', 'MRE11A', 'SCN1A', 'SDHA', 'KCNQ1', 'LMNA', 'ANK2', 'AXIN2', 'LAMA2', 'KCNH2', 'CACNA1C', 'ASPM', 'SMAD4', 'PCNT', 'MEN1', 'FLNA', 'ABCA4', 'MECP2', 'BAP1', 'PKP2', 'FLCN', 'MYO7A', 'BLM', 'APOB', 'ATP7B', 'VHL', 'MET', 'WFS1', 'DYNC1H1', 'CDKN2A', 'CNTNAP2', 'KMT2D', 'RELN', 'FANCC', 'VPS13B', 'NOTCH1',
                      'FH', 'ALMS1', 'DSG2', 'PKD1', 'SPTAN1', 'ZNF469', 'ACTN2', 'KCNT1', 'FLNC', 'CHD7', 'COL5A2', 'MYLK', 'CACNA1A', 'CAPN3', 'CACNA1S', 'RBM20', 'SCN4A']:
        symbol = 'Other'

    feature = feature[0:2]
    if feature not in ['NM', 'XM']:
        feature = 'Other'

    if exon not in ['16/16', '11/27', '4/10', '3/3', '2/2', '10/24', '326/363', '4/13', '4/4', '4/11', '8/8', '1/10', '9/10', '11/11', '1/1', '10/16', '11/15', '13/16', '10/10', '5/13', '12/16', '8/10', '3/16', '2/3', '1/3', '5/10', '6/10', '2/10', '2/16', '15/16', '7/10', '2/4', '14/16', '11/16', '9/9', '7/16', '28/28', '8/11', '9/16', '15/15', '5/16', '6/11', '10/27', '7/7', '4/16', '20/20', '6/6', '5/11', '3/10', '13/13', '7/11', '2/9', '33/33', '1/2', '4/18', '8/15', '5/6', '6/16', '12/13', '3/6', '1/9',
                    '12/15', '12/12', '7/8', '3/9', '5/5', '1/16', '358/363', '2/6', '1/4', '27/27', '2/5', '10/11', '3/4', '8/16', '10/13', '24/24', '7/9', '2/11', '5/8', '1/11', '1/6', '7/13', '4/7', '3/11', '23/24', '5/9', '4/9', '9/11', '5/15', '8/9', '32/33', '6/9', '9/13', '4/8', '1/5', '7/15', '14/15', '13/15', '4/6', '3/5', '2/12', '6/8', '2/13', '17/17', '1/8', '4/5', '9/20', '2/7', '11/13', '2/8', '3/7', '2/15', '3/22', '18/28', '8/13', '1/7', '5/12', '6/15', '6/13', '7/20', '5/7', '8/23', '18/27', '11/20', '15/20']:
        exon = 'Other'

    if amino_acids not in ['A', 'L', 'P', 'S', 'T', 'R/Q', 'R/H', 'G', 'R/C', 'A/T', 'D', 'R',
                           'P/L', 'E/K', 'V', 'A/V', 'R/W', 'V/I', 'N', 'I/V', 'R/*', 'V/M', 'I',
                           'Y', 'E', 'N/S', 'D/N', 'G/R', 'T/M', 'P/S', 'G/S', 'Q', 'H', 'I/T',
                           'T/I', 'F', 'T/A', 'K', 'S/L', 'M/V', 'K/R', 'Q/*', 'Y/C', 'L/F', 'L/V',
                           'Q/R', 'V/A', 'C', 'K/E', 'V/L', 'D/G', 'G/D', 'H/R', 'L/P', 'S/N',
                           'E/D', 'D/E', 'M/T', 'R/G', 'E/G', 'S/F', 'A/S', 'T/S', 'F/L', 'G/E',
                           'P/R', 'C/Y', 'Q/H', 'S/G', 'G/V', 'M/I', 'P/A', 'N/D', 'S/C', 'H/Y',
                           'E/Q', 'S/R', 'N/K', 'K/N', 'I/M', 'S/T', 'A/G', 'P/T', 'R/K', 'W/*',
                           'S/P', 'G/A', 'Y/*', 'Q/E', 'D/Y', 'C/R', 'R/L', 'R/S', 'D/H', '-/X',
                           'A/P', 'E/*', 'Y/H', 'D/V', 'R/P']:
        amino_acids = 'Other'

    CHROM_encoded = encoder_1.transform(np.array([chrom]))
    REF_encoded = encoder_2.transform(np.array([ref]))
    ALT_encoded = encoder_3.transform(np.array([alt]))
    CLNVC_encoded = encoder_4.transform(np.array([clnvc]))
    Allele_encoded = encoder_5.transform(np.array([allele]))
    IMPACT_encoded = encoder_6.transform(np.array([impact]))
    SYMBOL_encoded = encoder_7.transform(np.array([symbol]))
    Feature_encoded = encoder_8.transform(np.array([feature]))
    EXON_encoded = encoder_9.transform(np.array([exon]))
    Amino_acids_encoded = encoder_10.transform(np.array([amino_acids]))   

    t = [
    [
        float(CHROM_encoded), float(REF_encoded), float(ALT_encoded), float(af_esp), float(af_exac), float(af_tgp), float(CLNVC_encoded),
        float(clnvi), float(origin), float(Allele_encoded), float(IMPACT_encoded), 
        float(SYMBOL_encoded), float(Feature_encoded), float(EXON_encoded), float(intron), float(Amino_acids_encoded), float(strand), 
        float(bam_edit), float(sift), float(polyphen), float(loftool), float(cadd_phred), float(cadd_raw), 
        float(blosum62)
    ]
    ]

    t_scaled = scaler.transform(t)

    output = model.predict(t_scaled)
    print(output)

    return render_template('index.html', y="The prodicted answer is "+str(output[0]))


if __name__ == '__main__':
    app.run(debug=True)