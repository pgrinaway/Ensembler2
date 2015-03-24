__author__ = 'Patrick B. Grinaway'


import PDB
from MSMSeed import MSMSeed

def _read_local_pdb_repository(local_repo, pdb_code_input):
    """
    Utility function to read from the local PDB repository (this is considerably
    more reliable than trying to download relevant PDBs each time)

    Arguments
    ---------
    local_repo : String
        The location of the local repository
    pdb_code_input : String
        The PDB code of the structure to retrieve

    Returns
    -------
    A string containing the contents of the requested
    PDB file

    """
    import os
    import gzip
    pdb_code = pdb_code_input.lower()
    pdb_filename = 'pdb%s.ent.gz' % pdb_code
    pdb_path = '%s/%s' %(pdb_code[1:3], pdb_filename)
    filepath = os.path.join(local_repo, pdb_path)
    return "".join(gzip.GzipFile(filename=filepath).readlines())

def blast_pdb_local(fasta_string, num_hits=1000):
    """
    A utility function to generate the initial data for Ensembler from a local
    BLAST database as well as local PDB database, specified in the $PDB_HOME
    directory

    Arguments
    ---------
    fasta_string : String
        A FASTA-format string containing the sequence to search for
    num_hits : int, optional (default 1000)
        The maximum number of hits from the BLAST search

    Returns
    -------
    msmseeds : list of MSMSeed objects
        A list of MSMSeed objects that can be used in a distributed computing framework such as Apache Spark

    """
    import subprocess
    import os
    import shlex
    import StringIO
    import simtk.openmm.app as app
    blast_data = os.getenv("DATA_HOME")
    blast_query = 'blastp -db %s/pdbaa -max_target_seqs %d -outfmt' % (blast_data, num_hits)
    out_fmt = '7'
    blast_cmd = shlex.split(blast_query)
    blast_cmd.append(out_fmt)
    p = subprocess.Popen(blast_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    blast_aln, error = p.communicate(input=fasta_string)
    msmseeds = []
    local_pdb_repo = os.getenv("PDB_HOME")
    for i, result in enumerate(blast_aln.splitlines()):
        if result[0]!="#":
            res_data = result.split("\t")
            e_value = float(res_data[2])
            template_chain_code = "_".join(res_data[1].split("|")[3:])
            raw_template_pdb = _read_local_pdb_repository(local_pdb_repo, template_chain_code.split("_")[0])
            template_fasta, pdb_resnums = _retrieve_fasta(template_chain_code)
            template_pdb = StringIO.StringIO()
            raw_template_pdbio = StringIO.StringIO(raw_template_pdb)
            raw_template_pdbio.seek(0)
            end_resnums = PDB.extract_residues_by_resnum(template_pdb, raw_template_pdbio, pdb_resnums, template_chain_code.split("_")[1])
            template_pdb.seek(0)
            if template_pdb.len == 0:
                continue
            template_pdbfile = app.PDBFile(template_pdb)
            msmseeds.append(MSMSeed(fasta_string, template_fasta, template_pdbfile, e_value, model_id=i))
    return msmseeds


def retrieve_sifts(pdb_id):
    """Retrieves a SIFTS .xml file, given a PDB ID. Works by modifying the PDBe download URL.
    Also removes annoying namespace stuff.

    Arguments
    ---------
    pdb_id : String
        The ID of the PDB whose SIFTS should be retrieved

    Returns
    -------
    sifts_page_processed : String
        SIFTS xml with annoying namespace stuff deleted.
    """
    import re, gzip, StringIO, urllib2
    sifts_download_base_url='ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'
    url = sifts_download_base_url + pdb_id.lower() + '.xml.gz'
    try:
        response = urllib2.urlopen(url)
        sifts_page = response.read(100000000) # Max 100MB
    except:
        return ""
    # Decompress string
    sifts_page = gzip.GzipFile(fileobj=StringIO.StringIO(sifts_page)).read()

    # Removing all attribs from the entry tag, and the rdf tag and contents
    sifts_page_processed = ''
    skip_rdf_tag_flag = False
    for line in sifts_page.splitlines():
        if line[0:6] == '<entry':
            sifts_page_processed += '<entry>' + '\n'
        elif line[0:7] == '  <rdf:':
            skip_rdf_tag_flag = True
            pass
        elif line[0:8] == '  </rdf:':
            skip_rdf_tag_flag = False
            pass
        else:
            if skip_rdf_tag_flag:
                continue
            sifts_page_processed += line + '\n'
    return sifts_page_processed


def _retrieve_fasta(pdb_code_input):
    """
    This is a utility function to retrieve the FASTA of a PDB

    Arguments
    ---------
    pdb_code_input : String
        A string of format PDBCODE_CHAINCODE to select which PDB and chain should be retrieved

    Returns
    -------
    A FASTA-formatted string of the PDB chain

    The number of residues retrieved
    """
    from lxml import etree
    import StringIO
    pdb_code, chain_code = pdb_code_input.split("_")
    sifts = ""
    while sifts == "":
        sifts = retrieve_sifts(pdb_code)
    parser = etree.XMLParser(huge_tree=True)
    siftsXML = etree.parse(StringIO.StringIO(sifts), parser).getroot()
    observed_residues = siftsXML.xpath('entity/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"][not(../residueDetail[contains(text(),"Not_Observed")])][../crossRefDb[@dbSource="UniProt"]][not(../residueDetail[contains(text(),"modified")])][not(../residueDetail[contains(text(),"Conflict")])][not(../residueDetail[contains(text(),"mutation")])]' % chain_code)
    pdb_seq = ''.join([residue.find('../crossRefDb[@dbSource="UniProt"]').get('dbResName') for residue in observed_residues])
    template_PDBresnums = [residue.get('dbResNum') for residue in observed_residues]
    return "\n".join(['>' + pdb_code_input, pdb_seq]), template_PDBresnums



def _extract_seq(pdb):
    """
    Utility function to extract sequence. Deprecated.
    """
    import mdtraj as md
    mdtop = md.Topology.from_openmm(pdb.topology)
    resilist = []
    for residue in mdtop.residues:
        resilist.append(residue)
    return resilist












