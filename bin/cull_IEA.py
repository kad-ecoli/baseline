#!/usr/bin/env python
docstring='''cull_IEA.py go-basic.obo target.map goa_uniprot_all.gaf  goa_uniprot_all.clean.gaf goa_uniprot_all.is_a
    remove proteins and terms unrelated to CAFA dataset from goa_uniprot_all.gaf.
    save cleaned gaf file to goa_uniprot_all.is_a
'''
import sys
import obo2csv # parsing GO hierachy

def read_map_file(map_file):
    accession_list=[]
    fp=open(map_file,'r')
    for line in fp.read().splitlines():
        target,accession=line.split()
        accession_list.append(accession)
    fp.close()
    return accession_list

def cull_IEA(GOterm_list,accession_list,gaf_file,out_file):
    txt=''
    fp=open(gaf_file,'r')
    for line in fp:
        items=line.split('\t')
        accession=items[1]
        qualifier=items[3]
        GOterm=items[4]
        evidence=items[6]
        if  accession in accession_list and \
            not qualifier.startswith("NOT") and \
            GOterm in GOterm_list and evidence!="ND":
            txt+=line
    fp.close()
    fp=open(out_file,'w')
    fp.write(txt)
    fp.close()
    return

def isa_IEA(obo_dict,out_file,isa_file):
    annotation_dict=dict()
    fp=open(out_file,'r')
    for line in fp:
        items=line.split('\t')
        accession=items[1]
        GOterm=items[4]
        evidence=items[6]
        Aspect=items[8]
        key='\t'.join([accession,Aspect,GOterm])
        if not key in annotation_dict:
            annotation_dict[key]=[]
        elif evidence in annotation_dict[key]:
            continue
        annotation_dict[key].append(evidence)
    fp.close()
    for key,evidence_list in annotation_dict.items():
        accession,Aspect,GOterm=key.split('\t')
        for parent_GO in obo_dict.is_a(Term_id=GOterm, direct=False,
            name=True, number=False).split('\t'):
            GO_ID,name=parent_GO.split(" ! ")
            key='\t'.join([accession,Aspect,GO_ID])
            if not key in annotation_dict:
                annotation_dict[key]=evidence_list
            else:
                annotation_dict[key]+=evidence_list
    for key in annotation_dict:
        annotation_dict[key]=sorted(set(annotation_dict[key]))
    txt=''
    for key in sorted(annotation_dict.keys()):
        txt+=key+'\t'+','.join(annotation_dict[key])+'\n'
    fp=open(isa_file,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=6:
        sys.stderr.write(docstring)
        exit()

    obo_file=sys.argv[1]
    map_file=sys.argv[2]
    gaf_file=sys.argv[3]
    out_file=sys.argv[4]
    isa_file=sys.argv[5]

    obo_dict=obo2csv.parse_obo_file(obo_file)
    accession_list=read_map_file(map_file)
    GOterm_list=obo_dict['F']['Term'].keys()+ \
                obo_dict['P']['Term'].keys()+ \
                obo_dict['C']['Term'].keys()
    cull_IEA(GOterm_list,accession_list,gaf_file,out_file)
    isa_IEA(obo_dict,out_file,isa_file)
