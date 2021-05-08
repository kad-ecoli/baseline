#include "basic_fun.h"
#include "BLOSUM.h"

using namespace std;

void print_extra_help()
{
    cout <<
"Additional options:\n"
"    -atom    4-character atom name used to represent a residue.\n"
"             Default is \" C3'\" for RNA/DNA and \" CA \" for proteins\n"
"             (note the spaces before and after CA).\n"
"\n"
"    -mol     Molecule type: RNA or protein\n"
"             Default is detect molecule type automatically\n"
"\n"
"    -ter     Strings to mark the end of a chain\n"
"             3: TER, ENDMDL, END or different chain ID\n"
"             2: ENDMDL, END, or different chain ID\n"
"             1: ENDMDL or END\n"
"             0: (default) end of file\n"
"\n"
"             For FASTA intput (-infmt 4), -ter 0 means read all\n"
"             sequences; -ter >=1 means read the first sequence only."
"\n"
"    -split   Whether to split PDB file into multiple chains\n"
"             0: treat the whole structure as one single chain\n"
"             1: treat each MODEL as a separate chain (-ter should be 0)\n"
"             2: (default) treat each chain as a seperate chain (-ter should be <=1)\n"
"\n"
"             For FASTA intput, -split 0 means concatenate all sequences into\n"
"             one read all sequence; -split >=1 means each sequence is an\n"
"             individual entry."
"\n"
"    -het     Whether to align residues marked as 'HETATM' in addition to 'ATOM  '\n"
"             0: (default) only align 'ATOM  ' residues\n"
"             1: align both 'ATOM  ' and 'HETATM' residues\n"
    <<endl;
}

void print_help(bool h_opt=false)
{
    cout <<
"Use BLOSUM62 to get the raw score and bit score of sequences aligning to themselves.\n"
"\n"
"Usage: SelfScore seq.fasta\n"
"\n"
"Options:\n"
"    -h       Print the full help message\n"
"\n"
"    -infmt   Input format\n"
"            -1: automatically detect PDB or PDBx/mmCIF format\n"
"             0: PDB format\n"
"             3: PDBx/mmCIF format\n"
"             4: (default) FASTA format sequence\n"
    <<endl;

    if (h_opt) print_extra_help();

    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    if (argc < 2) print_help();

    /**********************/
    /*    get argument    */
    /**********************/
    string xname     ="";
    string yname     ="";
    bool   h_opt     =false; // print full help message
    int    infmt_opt =4;     // FASTA sequence
    int    ter_opt   =0;     // end of file
    int    split_opt =2;     // split chain
    int    het_opt=0;        // do not read HETATM residues
    string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA

    for(int i = 1; i < argc; i++)
    {
        if ( !strcmp(argv[i],"-h") )
        {
            h_opt = true;
        }
        else if ( !strcmp(argv[i],"-infmt") && i < (argc-1) )
        {
            infmt_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-ter") && i < (argc-1) )
        {
            ter_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-split") && i < (argc-1) )
        {
            split_opt=atoi(argv[i + 1]); i++;
        }
        else if ( !strcmp(argv[i],"-atom") && i < (argc-1) )
        {
            atom_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-mol") && i < (argc-1) )
        {
            mol_opt=argv[i + 1]; i++;
        }
        else if ( !strcmp(argv[i],"-het") && i < (argc-1) )
        {
            het_opt=atoi(argv[i + 1]); i++;
        }
        else if (xname.size() == 0) xname=argv[i];
        else PrintErrorAndQuit(string("ERROR! Undefined option ")+argv[i]);
    }

    if(xname.size()==0)
    {
        if (h_opt) print_help(h_opt);
        if (xname.size()==0)
            PrintErrorAndQuit("Please provide input sequences");
    }
    if (atom_opt.size()!=4)
        PrintErrorAndQuit("ERROR! Atom name must have 4 characters, including space.");
    if (mol_opt!="auto" && mol_opt!="protein" && mol_opt!="RNA")
        PrintErrorAndQuit("ERROR! Molecule type must be either RNA or protein.");
    else if (mol_opt=="protein" && atom_opt=="auto")
        atom_opt=" CA ";
    else if (mol_opt=="RNA" && atom_opt=="auto")
        atom_opt=" C3'";

    if (split_opt==1 && ter_opt!=0)
        PrintErrorAndQuit("-split 1 should be used with -ter 0");
    else if (split_opt==2 && ter_opt!=0 && ter_opt!=1)
        PrintErrorAndQuit("-split 2 should be used with -ter 0 or 1");
    if (split_opt<0 || split_opt>2)
        PrintErrorAndQuit("-split can only be 0, 1 or 2");

    /* parse file list */
    cout<<"#target\tlength\tscore\tbitscore"<<endl;

    /* declare previously global variables */
    vector<vector<string> >PDB_lines1; // text of chain1
    vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    vector<string> chainID_list1;      // list of chainID1
    int  i,j;                // file index
    int  chain_i,chain_j;    // chain index
    int  xlen, ylen;         // chain length
    int  xchainnum,ychainnum;// number of chains in a PDB file
    char *seqx, *seqy;       // for the protein sequence 
    int  l;                  // residue index

    /* Karlin-Altschul parameters */
    /* hsp->score = BLAST_Nint(((double) hsp->score) / scoreDivisor);
     * hsp->bit_score = (hsp->score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;
     * logK=log(kappa)
     * NCBIMATH_LN2=log(2)
     */
    /* ungapped parameters for BLOSUM62 */
    //const double lambda=0.3176;
    //const double kappa=0.134;
    /* gapped parameters for BLOSUM62 */
    //const double lambda=0.251;
    //const double kappa=0.031;
    /* current blastp implement */
    const double lambda=0.267;
    const double kappa=0.041;

    /* loop over file names */
    if (infmt_opt>=4) xchainnum=get_FASTA_lines(xname, PDB_lines1, 
            chainID_list1, mol_vec1, ter_opt, split_opt);
    else xchainnum=get_PDB_lines(xname, PDB_lines1, chainID_list1,
            mol_vec1, ter_opt, infmt_opt, atom_opt, split_opt, het_opt);
    if (!xchainnum)
    {
        cerr<<"Warning! Cannot parse file: "<<xname
            <<". Chain number 0."<<endl;
        return 1;
    }
    for (chain_i=0;chain_i<xchainnum;chain_i++)
    {
        if (infmt_opt>=4) xlen=PDB_lines1[chain_i][0].size();
        else xlen=PDB_lines1[chain_i].size();
        if (mol_opt=="RNA") mol_vec1[chain_i]=1;
        else if (mol_opt=="protein") mol_vec1[chain_i]=-1;
        if (!xlen)
        {
            cerr<<"Warning! Cannot parse file: "<<xname
                <<". Chain length 0."<<endl;
            continue;
        }
        seqx = new char[xlen + 1];
        if (infmt_opt>=4) strcpy(seqx,PDB_lines1[chain_i][0].c_str());
        else for (l=0;l<xlen;l++)
            seqx[l]=AAmap(PDB_lines1[chain_i][l].substr(17,3));
        seqx[xlen]=0;

        long aln_score=0;
        for (l=0;l<xlen;l++)
            aln_score+=BLOSUM[seqx[l]][seqx[l]];
        double bit_score=(lambda*aln_score-log(kappa))/log(2);
        cout<<chainID_list1[chain_i]<<'\t'<<xlen<<'\t'<<aln_score<<'\t'
            <<setiosflags(ios::fixed)<<setprecision(1)<<bit_score<<endl;
        
        PDB_lines1[chain_i].clear();
        delete [] seqx;
    } // chain_i
    xname.clear();
    PDB_lines1.clear();
    chainID_list1.clear();
    mol_vec1.clear();
    return 0;
}
