import sys, os, random, copy

import gffutils

def mod_transcript(tr, i):
    i = str(i)
    new_tr_name = tr.attributes['transcript_name'][0] + "_" + i
    new_tr_id = tr.attributes['transcript_id'][0] + "_" + i
    try:
        new_tr_hav = tr.attributes['havana_transcript'][0] + "_" + i
    except KeyError:
        new_tr_hav = ""
    tr1 = copy.deepcopy(tr)
    tr1.attributes['transcript_name'][0] = new_tr_name
    tr1.attributes['transcript_id'][0] = new_tr_id
    if new_tr_hav != "":
        tr1.attributes['havana_transcript'][0] = new_tr_hav
    return tr1

def mod_exon(ex, i):
    i = str(i)
    ex_tr_id = ex.attributes['transcript_id'][0] + "_" + i
    ex_tr_name = ex.attributes['transcript_name'][0] + "_" + i
    try:
        ex_tr_hav = ex.attributes['havana_transcript'][0] + "_" + i
    except KeyError:
        ex_tr_hav = ""
    ex1 = copy.deepcopy(ex)
    ex1.attributes['transcript_id'][0] = ex_tr_id
    ex1.attributes['transcript_name'][0] = ex_tr_name
    if ex_tr_hav != "":
        ex1.attributes['havana_transcript'][0] = ex_tr_hav
    return ex1

def compiting(ex, i, start, end):
    i = str(i)
    ex_s = ex.start
    ex_e = ex.end
    if start and end:
        if random.randint(0,100) < 50:
            ex_s += 12
        else:
            ex_e -= 12
    else:
        if start:
            ex_s += 12
        if end:
            ex_e -= 12
    ex_tr_id = ex.attributes['transcript_id'][0] + "_" + i
    ex_tr_name = ex.attributes['transcript_name'][0] + "_" + i
    ex_id = ex.attributes['exon_id'][0] + "_" + i
    try:
        ex_tr_hav = ex.attributes['havana_transcript'][0] + "_" + i
    except KeyError:
        ex_tr_hav = ""
    ex1 = copy.deepcopy(ex)
    ex1.start = ex_s
    ex1.end = ex_e
    ex1.attributes['transcript_id'][0] = ex_tr_id
    ex1.attributes['transcript_name'][0] = ex_tr_name
    ex1.attributes['exon_id'][0] = ex_id
    if ex_tr_hav != "":
        ex1.attributes['havana_transcript'][0] = ex_tr_hav
    return ex1

def getExID(x,y):
    return str(x) + "_" + str(y)

def main():
    gtf_path = sys.argv[1]
    P = int(sys.argv[2])
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtf_path,
                                 dbfn="{}.db".format(gtf_path),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
    out = open("{}.new.gtf".format(gtf_path), "w")
    log = open("{}.log".format(gtf_path), "w")
    transcripts = []
    forbidden_exons = []
    for gene in gtf.features_of_type('gene'):
        for tr in gtf.children(gene, featuretype='transcript', order_by='start'):
            transcript = []
            for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                transcript.append(getExID(ex.start, ex.end))
            transcripts.append(transcript)

    for gene in gtf.features_of_type('gene'):
        log.write('Gene {}\n'.format(gene.id))
        out.write(str(gene) + "\n")
        tr_index = -1
        for tr in gtf.children(gene, featuretype='transcript', order_by='start'):
            tr_index += 1
            log.write('Transcript {}\n'.format(tr.id))
            #out.write(str(tr) + "\n")
            exons = []
            for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                exons.append(ex)

            if len(exons) < 3:
                out.write(str(tr) + "\n")
                for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                    out.write(str(ex) + "\n")
            else:
                mod_flag = True
                mod_number = 1
                # --------------------------------------------------------------------
                # Exon Skipping
                if random.randint(0,100) < P:
                    mod_flag = False
                    tr1 = mod_transcript(tr, mod_number)
                    log.write('Transcript {}_{}\n'.format(tr.id, mod_number))
                    out.write(str(tr1) + "\n")
                    ex_to_skip_i = random.randint(0,len(exons)-1)
                    log.write('- ES: Exon {} ({})\n'.format(exons[ex_to_skip_i].id, ex_to_skip_i+1))
                    i = 0
                    for ex in exons:
                        if i != ex_to_skip_i:
                            ex1 = mod_exon(ex, mod_number)
                            if i == ex_to_skip_i-1:
                                if random.randint(0,100) < P:
                                    ex1 = compiting(ex, mod_number, False, True)
                                    log.write('- Comp: Exon {} ({})\n'.format(exons[i].id, i+1))
                            if i == ex_to_skip_i+1:
                                if random.randint(0,100) < P:
                                    log.write('- Comp: Exon {} ({})\n'.format(exons[i].id, i+1))
                                    ex1 = compiting(ex, mod_number, True, False)
                            out.write(str(ex1) + "\n")
                        i+=1
                    mod_number += 1
                # --------------------------------------------------------------------
                # Mutually Exclusive Exons
                if random.randint(0,100) < P:
                    mod_flag = False
                    flag = True
                    ids = list(range(0,len(exons)-1))
                    while flag and ids != []:
                        flag = False
                        ex_i1 = random.choice(ids)
                        ids.remove(ex_i1)
                        ex_i2 = ex_i1+1
                        i=-1
                        for transcript in transcripts:
                            i+=1
                            if i!=tr_index:
                                ex1 = exons[ex_i1]
                                ex2 = exons[ex_i2]
                                getExID(ex.start, ex.end)
                                if getExID(ex1.start, ex1.end) in transcript and getExID(ex2.start, ex2.end) in transcript:
                                    flag = True
                                    break
                    if not(flag):
                        tr1 = mod_transcript(tr, mod_number)
                        log.write('Transcript {}_{}\n'.format(tr.id, mod_number))
                        out.write(str(tr1) + "\n")
                        log.write('- MEE: Exons {} ({}) {} ({})\n'.format(exons[ex_i1].id, ex_i1+1, exons[ex_i2].id, ex_i2+1))
                        i = 0
                        for ex in exons:
                            if i != ex_i1:
                                ex1 = mod_exon(ex, mod_number)
                                if i == ex_i1-1:
                                    if random.randint(0,100) < P:
                                        ex1 = compiting(ex, mod_number, False, True)
                                        log.write('- Comp: Exon {} ({})\n'.format(exons[i].id, i+1))
                                if i == ex_i1+1:
                                    if random.randint(0,100) < P:
                                        log.write('- Comp: Exon {} ({})\n'.format(exons[i].id, i+1))
                                        ex1 = compiting(ex, mod_number, True, False)
                                out.write(str(ex1) + "\n")
                            i+=1
                        mod_number += 1
                        tr1 = mod_transcript(tr, mod_number)
                        out.write(str(tr1) + "\n")
                        i = 0
                        for ex in exons:
                            if i != ex_i2:
                                ex1 = mod_exon(ex, mod_number)
                                if i == ex_i2-1:
                                    if random.randint(0,100) < P:
                                        ex1 = compiting(ex, mod_number, False, True)
                                        log.write('- Comp: Exon {} ({})\n'.format(exons[i].id, i+1))
                                if i == ex_i2+1:
                                    if random.randint(0,100) < P:
                                        log.write('- Comp: Exon {} ({})\n'.format(exons[i].id, i+1))
                                        ex1 = compiting(ex, mod_number, True, False)
                                out.write(str(ex1) + "\n")
                            i+=1
                        mod_number += 1
                # --------------------------------------------------------------------
                # Compiting
                if random.randint(0,100) < P:
                    mod_flag = False
                    tr1 = mod_transcript(tr, mod_number)
                    log.write('Transcript {}_{}\n'.format(tr.id, mod_number))
                    out.write(str(tr1) + "\n")
                    i = 0
                    ex_to_comp_i = random.randint(0,len(exons)-1)
                    log.write('- Comp: Exon {} ({})\n'.format(exons[ex_to_comp_i].id, ex_to_comp_i+1))
                    for ex in exons:
                        if i == ex_to_comp_i:
                            ex1 = compiting(ex, mod_number, True, True)
                            out.write(str(ex1) + "\n")
                        else:
                            ex1 = mod_exon(ex, mod_number)
                            out.write(str(ex1) + "\n")
                        i+=1
                if mod_flag:
                    out.write(str(tr) + "\n")
                    exons = []
                    for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                        exons.append(ex)
                        out.write(str(ex) + "\n")
    out.close()
    log.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Given a GTF (maybe also GFF), this script produces a new GTF",
              "containing new transcripts produced by some alternative",
              "splicing events (exon skipping, mutually exclusive exons,",
              "competing). The events are added with probability P \\in [0,100].",
              "",
              "\tUsage: python3 modGTFs.py /path/to/GTF P",
              "",
              sep="\n")
    else:
        main()
