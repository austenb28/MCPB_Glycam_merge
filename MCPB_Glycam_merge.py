#!/usr/bin/env python3

import argparse
from array import array

def file_to_str(filename):
	with(open(filename, 'r')) as myfile:
		buff=myfile.read()
	return buff

def skipline(buff,index):
	while(buff[index] != '\n'):
		index += 1
	index += 1
	return index

def get_MCPB_pdb_str(leap_str):
	identifier = "mol = loadpdb "
	j = leap_str.index(identifier) + len(identifier)
	k = leap_str.index('\n',j)
	return file_to_str(leap_str[j:k])

def ASNS_to_NLNS(ASN_pdb_str):
	j = 0
	N = len(ASN_pdb_str)
	NLN_pdb_list = []
	patterns = []
	while j < N:
		if(ASN_pdb_str[j+17:j+20] == 'ASN'):
			patterns.append(ASN_pdb_str[j+17:j+26])
		else:
			break
		j = skipline(ASN_pdb_str,j)
	j = 0
	while j < N:
		k = skipline(ASN_pdb_str,j)
		if(ASN_pdb_str[j+17:j+26] in patterns):
			line = (
				ASN_pdb_str[j:j+17] + 
				'NLN' + 
				ASN_pdb_str[j+20:k]
			)
		else:
			line = ASN_pdb_str[j:k]
		NLN_pdb_list.append(line)
		j = k
	return ''.join(NLN_pdb_list)

def get_natoms_nres(pdb_str):
	j = 0
	N = len(pdb_str)
	while j+11 < N:
		if (pdb_str[j:j+4] == 'ATOM' or
			pdb_str[j:j+6] == 'HETATM'):
			natoms = int(pdb_str[j+6:j+11])
			nres = int(pdb_str[j+22:j+26])
		j = skipline(pdb_str,j)
	return [natoms,nres]

def reindex_strs(leap_str,MCPB_pdb_str,GLYCAM_pdb_str):
	[natoms_GLYCAM,nres_GLYCAM] = get_natoms_nres(GLYCAM_pdb_str)
	[natoms_MCPB,nres_MCPB] = get_natoms_nres(MCPB_pdb_str)
	idmap_GLYCAM = array('I',[0]*(natoms_GLYCAM+1))
	idmap_MCPB = array('I',[0]*(natoms_MCPB+1))
	resmap_GLYCAM = array('I',[0]*(nres_GLYCAM+1))
	resmap_MCPB = array('I',[0]*(nres_MCPB+1))
	final_pdb_lst = []
	j = GLYCAM_pdb_str.index('ATOM')
	N_j = len(GLYCAM_pdb_str)
	k = MCPB_pdb_str.index('ATOM')
	N_k = len(MCPB_pdb_str)
	atomid = 1
	resid = 1
	line = ''
	j2 = skipline(GLYCAM_pdb_str,j)
	k2 = skipline(MCPB_pdb_str,k)
	prev_resid_GLYCAM = 0
	TER_count = int(0)
	while (j < N_j and k < N_k):
		# print(GLYCAM_pdb_str[j:j2])
		# print(MCPB_pdb_str[k:k2])
		if (GLYCAM_pdb_str[j:j+3] == 'TER' and 
			MCPB_pdb_str[k:k+3] == 'TER'):
			final_pdb_lst.append(GLYCAM_pdb_str[j:j2])
			j = j2
			k = k2
			TER_count += 1
			j2 = skipline(GLYCAM_pdb_str,j)
			k2 = skipline(MCPB_pdb_str,k)
		elif (GLYCAM_pdb_str[j:j+6] == 'HETATM'):
			while(GLYCAM_pdb_str[j:j+6] == 'HETATM' or
				GLYCAM_pdb_str[j:j+3] == 'TER'):
				if GLYCAM_pdb_str[j:j+3] == 'TER':
					final_pdb_lst.append(GLYCAM_pdb_str[j:j2])
					j = j2
					j2 = skipline(GLYCAM_pdb_str,j)
				else:
					resid_GLYCAM = GLYCAM_pdb_str[j+22:j+26]
					if (prev_resid_GLYCAM != 0 and
						prev_resid_GLYCAM != resid_GLYCAM):
						resid += 1
					final_pdb_lst.append(GLYCAM_pdb_str[j:j+6] +
						'{:>5d}'.format(atomid) +
						GLYCAM_pdb_str[j+11:j+22] +
						'{:>4d}'.format(resid) +
						GLYCAM_pdb_str[j+26:j2]
					)
					resmap_GLYCAM[int(resid_GLYCAM)] = resid
					ind = int(GLYCAM_pdb_str[j+6:j+11])
					idmap_GLYCAM[ind] = atomid
					j = j2
					j2 = skipline(GLYCAM_pdb_str,j)
					atomid += 1
					prev_resid_GLYCAM = resid_GLYCAM
			resid += 1
		elif(MCPB_pdb_str[k:k+3] == 'TER'):
			TER_count += 1
			# print(MCPB_pdb_str[k:k2])
			final_pdb_lst.append(MCPB_pdb_str[k:k2])
			k = k2
			k2 = skipline(MCPB_pdb_str,k)
		elif(GLYCAM_pdb_str[j:j+3] == 'TER'):
			j = j2
			j2 = skipline(GLYCAM_pdb_str,j)
		elif(MCPB_pdb_str[k:k+3] != 'END' and
			(MCPB_pdb_str[k:k+4] == 'ATOM' or 
			MCPB_pdb_str[k:k+6] == 'HETATM')
		):
			if(GLYCAM_pdb_str[j+17:j+20] == MCPB_pdb_str[k+17:k+20] or
				GLYCAM_pdb_str[j+17:j+20] == 'NLN'):
				resid_GLYCAM = GLYCAM_pdb_str[j+22:j+26]
				cur_resid = resid_GLYCAM
				resmap_GLYCAM[int(cur_resid)] = resid
				while (cur_resid == resid_GLYCAM and 
					GLYCAM_pdb_str[j:j+3] != 'TER'
				):
					if(GLYCAM_pdb_str[j+17:j+20] == 'CYX' and
						GLYCAM_pdb_str[j+13:j+15] == 'SG'
					):
						SG_id = atomid
					final_pdb_lst.append(GLYCAM_pdb_str[j:j+4] +
						'{:>7d}'.format(atomid) +
						GLYCAM_pdb_str[j+11:j+22] +
						'{:>4d}'.format(resid) +
						GLYCAM_pdb_str[j+26:j2]
					)
					ind = int(GLYCAM_pdb_str[j+6:j+11])
					idmap_GLYCAM[ind] = atomid
					j = j2
					j2 = skipline(GLYCAM_pdb_str,j)
					cur_resid = GLYCAM_pdb_str[j+22:j+26]
					atomid += 1
				resid_MCPB = MCPB_pdb_str[k+22:k+26]
				cur_resid = resid_MCPB
				resmap_MCPB[int(resid_MCPB)] = resid
				while (cur_resid == resid_MCPB and 
					MCPB_pdb_str[k:k+3] != 'TER'
				):
					if(MCPB_pdb_str[k+17:k+20] == 'CYX' and
						MCPB_pdb_str[k+13:k+15] == 'SG'
					):
						ind = int(MCPB_pdb_str[k+6:k+11]) - TER_count # doing this for pdb4amber, may be invalid if not using pdb4amber and using SG - SG CONECT
						idmap_MCPB[ind] = SG_id
					k = k2
					k2 = skipline(MCPB_pdb_str,k)
					cur_resid = MCPB_pdb_str[k+22:k+26]
				resid += 1
			else:
				# print(MCPB_pdb_str[k:k2])
				resid_MCPB = MCPB_pdb_str[k+22:k+26]
				cur_resid = resid_MCPB
				resmap_MCPB[int(resid_MCPB)] = resid
				while (cur_resid == resid_MCPB and 
					MCPB_pdb_str[k:k+3] != 'TER'
				):
					final_pdb_lst.append(MCPB_pdb_str[k:k+4] +
						'{:>7d}'.format(atomid) +
						MCPB_pdb_str[k+11:k+21] +
						' ' +
						'{:>4d}'.format(resid) +
						MCPB_pdb_str[k+26:k2]
					)
					ind = int(MCPB_pdb_str[k+6:k+11])
					k = k2
					k2 = skipline(MCPB_pdb_str,k)
					cur_resid = MCPB_pdb_str[k+22:k+26]
					atomid += 1
				resid += 1
		elif(MCPB_pdb_str[k:k+3] != 'END' and
			MCPB_pdb_str[k:k+6] == 'CONECT'
		):
			while (k < N_k):
				if (MCPB_pdb_str[k:k+3] == 'END'):
					break
				else:
					line = ['CONECT']
					kw = 5
					n = 1
					fs = k+6
					fe = fs + kw
					while n <= 2:
						field = MCPB_pdb_str[fs:fe]
						line.append('{:>5d}'.format(idmap_MCPB[int(field)]))
						fs += kw
						fe += kw
						n += 1
					line.append('\n')
					final_pdb_lst.append(''.join(line))
					k = k2
					k2 = skipline(MCPB_pdb_str,k)
		elif (GLYCAM_pdb_str[j:j+6] == 'CONECT'):
			# print(GLYCAM_pdb_str[j:j2])
			while (j < N_j):
				if j == N_j - 1:
					final_pdb_lst.append('\n')
					break
				elif GLYCAM_pdb_str[j:j+3] == 'END':
					final_pdb_lst.append(GLYCAM_pdb_str[j:j2])
					j = j2
				else:
					field = GLYCAM_pdb_str[j+6:j+11]
					line = ['CONECT']
					jw = 5
					n = 1
					while field != '     ':
						line.append('{:>5d}'.format(idmap_GLYCAM[int(field)]))
						n += 1
						field = GLYCAM_pdb_str[j+1+n*jw:j+1+(n+1)*jw]
					line.append(GLYCAM_pdb_str[j+1+n*jw:j2])
					final_pdb_lst.append(''.join(line))
					j = j2
					j2 = skipline(GLYCAM_pdb_str,j)
		elif (GLYCAM_pdb_str[j:j+3] == 'END' and
			MCPB_pdb_str[k:k+3] == 'END'
		):
			final_pdb_lst.append(GLYCAM_pdb_str[j:j2])
			break
	j = 0
	j2 = skipline(GLYCAM_pdb_str,j)
	link_list = []
	while (GLYCAM_pdb_str[j:j+4] == 'LINK'):
		res1 = int(GLYCAM_pdb_str[j+22:j+26])
		res2 = int(GLYCAM_pdb_str[j+52:j+56])
		res1 = resmap_GLYCAM[res1]
		res2 = resmap_GLYCAM[res2]
		link_list.append(
			GLYCAM_pdb_str[j:j+22] +
			'{:>4d}'.format(res1) + 
			GLYCAM_pdb_str[j+26:j+52] +
			'{:>4d}'.format(res2) +
			GLYCAM_pdb_str[j+56:j2]
		)
		j = j2
		j2 = skipline(GLYCAM_pdb_str,j)

	final_leapin_lst = []
	try:
		j = leap_str.index('source leaprc.gaff')
		j = skipline(leap_str,j)
		j2 = leap_str.index('mol = ',j)
		final_leapin_lst.append(
			leap_str[0:j] +
			'source leaprc.GLYCAM_06j-1\n' +
			leap_str[j:j2] +
			'mol = loadpdb merged.pdb\n'
		)
		j = leap_str.index('bond mol.',j2)
		j2 = skipline(leap_str,j)
		N = len(leap_str)
		while(leap_str[j:j+4] == 'bond'):
			k1 = leap_str.index('.',j)+1
			k2 = leap_str.index('.',k1)
			k3 = leap_str.index('.',k2+1)+1
			k4 = leap_str.index('.',k3)
			final_leapin_lst.append(
				leap_str[j:k1] +
				str(resmap_MCPB[int(leap_str[k1:k2])]) +
				leap_str[k2:k3] +
				str(resmap_MCPB[int(leap_str[k3:k4])]) +
				leap_str[k4:j2]
			)
			j = j2
			j2 = skipline(leap_str,j)
		final_leapin_lst.append(
			"saveamberparm mol merged.prmtop merged.inpcrd\n" +
			"quit\n"
		)
	except:
		print("\nWARNING: Ignoring MCPB leap input.  Either something's wrong, or you're doing something non-standard.\n")

	return [''.join(link_list+final_pdb_lst),''.join(final_leapin_lst)]

def str_to_file(buff,filename):
	with open(filename,'w') as myfile:
		myfile.write(buff)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('leapin', help= 'MCPB leap input file (.in)')
	parser.add_argument('glypdb', help='GLYCAM pdb file (.pdb)')
	args = parser.parse_args()
	leap_str = file_to_str(args.leapin)
	MCPB_pdb_str = get_MCPB_pdb_str(leap_str)
	GLYCAM_pdb_str = file_to_str(args.glypdb)
	GLYCAM_pdb_str = ASNS_to_NLNS(GLYCAM_pdb_str)
	[final_pdb_str,final_leapin_str] = reindex_strs(
		leap_str,
		MCPB_pdb_str,
		GLYCAM_pdb_str
	)
	str_to_file(final_pdb_str,'merged.pdb')
	str_to_file(final_leapin_str,'merged_leap.in')

main()