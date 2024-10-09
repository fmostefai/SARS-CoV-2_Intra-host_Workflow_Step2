from scipy.stats import binom
from pathlib import Path
from collections import namedtuple
import collections
import sys
import re


def parse_fieldBase(fields):
    ''' function to parse each field base into base nucleotide and number of counts for that nucleotide'''
    # input: field containing nuc and base sperated by ':'
    # output:  ReadCountBase namedtuple filled
    bases = [
        ReadCountBase('A', 
                      int(fields[0]) + int(fields[4]), 
                      int(fields[0]), 
                      int(fields[4]))
        
        ,ReadCountBase('T', 
                      int(fields[1]) + int(fields[5]), 
                      int(fields[1]), 
                      int(fields[5]))
        
        ,ReadCountBase('C', 
                      int(fields[2]) + int(fields[6]), 
                      int(fields[2]), 
                      int(fields[6]))
        
        ,ReadCountBase('G', 
                      int(fields[3]) + int(fields[7]), 
                      int(fields[3]), 
                      int(fields[7]))
    ]

    return bases

def compute_indel_depth(indel_field):
    if indel_field == 'NA':
        return 0, 0
    else:
        delimiters = "|"
        regex_pattern = '|'.join(map(re.escape, delimiters))
        indels = {}
        for i in re.split(regex_pattern, indel_field):
            indels[i.split(':')[1]] = int(i.split(':')[0])
        
        plus_strand_count = []
        minus_strand_count = []
        
        for key, value in indels.items():
            if key.isupper():
                plus_strand_count.append(value)
            elif key.islower:
                minus_strand_count.append(value)
        
        return sum(plus_strand_count), sum(minus_strand_count)

def parse_line(line):
    '''function to parse line seperated by tab'''
    # input: input_file line
    # ouput: ReadCount namedtuple filled
    fields = line.rstrip().split('\t')
    bases = parse_fieldBase(fields[3:11])

    insertion_plus_strand_count,  insertion_minus_strand_count = compute_indel_depth(fields[11])
    deletion_plus_strand_count,  deletion_minus_strand_count = compute_indel_depth(fields[12])
    insertion_depth = insertion_plus_strand_count + insertion_minus_strand_count
    deletion_depth = deletion_plus_strand_count + deletion_minus_strand_count
    
    bases_dict = {base.base: base for base in bases}
    bases_dict = collections.OrderedDict(sorted(bases_dict.items(), key=lambda b: b[1].count, reverse=True))
    
    depth = (
            bases_dict['A'].count + 
            bases_dict['T'].count + 
            bases_dict['C'].count +
            bases_dict['G'].count
            )

    total_plus_count = (
                        bases_dict['A'].num_plus_strand + 
                        bases_dict['T'].num_plus_strand + 
                        bases_dict['C'].num_plus_strand + 
                        bases_dict['G'].num_plus_strand
                        )

    total_minus_count = (
                        bases_dict['A'].num_minus_strand + 
                        bases_dict['T'].num_minus_strand + 
                        bases_dict['C'].num_minus_strand + 
                        bases_dict['G'].num_minus_strand
                        )

    return ReadCount(*fields[:3], depth, bases_dict, total_plus_count, total_minus_count, insertion_depth, deletion_depth)

def iter_readcounts(file):
    ''' Iterator to read one line at a time'''
    with file.open() as f:
        next(f)
        for line in f:
            yield parse_line(line)


def alternativeprob(k,n,p):
    #print(k, n, p)
    if n!=0 and k/n < p:
         return binom.cdf(k, n, p)
    else:
        return 1

def worstprobstrandbias(nam,tmc,nap,tpc,vaf):
    #p=(nam+nap)/(tmc+tpc)
    #p=min(0.05,p)
    p=vaf
    return min(alternativeprob(nam, tmc, p),
               alternativeprob(nap, tpc, p))

# Namedtuple for  parsing fields of read count files
ReadCountTuple = namedtuple('ReadCountTuple', [
                    'chr',
                    'position',
                    'reference_base',
                    'depth',
                    'bases',
                    'total_plus_count',
                    'total_minus_count',
                    'insertion_depth',
                    'deletion_depth'
])

# Name tuple for parsing each field
ReadCountBase = namedtuple('ReadCountBase', [
                    'base',
                    'count',
                    'num_plus_strand',
                    'num_minus_strand'
                    ])

class ReadCount(ReadCountTuple):
    def alt_base(self):
        alt_base_reads = [(k,v) for (k,v) in self.bases.items() if k != self.reference_base]
        return alt_base_reads[0][1]

    def ref_base_plus_strand(self):
        for k,v in self.bases.items():
            if k == self.reference_base:
                return v.num_plus_strand
    
    def ref_base_minus_strand(self):
        for k,v in self.bases.items():
            if k == self.reference_base:
                return v.num_minus_strand

    def vaf(self):
        if  (int(self.depth)) > 0:
            return self.alt_base().count / (int(self.depth))
        else: return 0

    def allele_count(self):
        return len([(k,v) for (k,v) in self.bases.items() if v.count != 0])

    def plus_ratio(self):
        if (self.depth) > 0:
            return self.total_plus_count / self.depth
        else: return None

    def minus_ratio(self):
        if (self.depth) > 0:
            return self.total_minus_count / self.depth
        else: return None

    def alt_plus_vaf(self):
        if (self.total_plus_count) > 0:
            return self.alt_base().num_plus_strand / self.total_plus_count
        else: return None

    def alt_minus_vaf(self):
        if (self.total_minus_count) > 0:
            return self.alt_base().num_minus_strand / self.total_minus_count
        else: return None

    def p_value(self):
        return worstprobstrandbias(self.alt_base().num_minus_strand,
                                   self.total_minus_count,
                                   self.alt_base().num_plus_strand,
                                   self.total_plus_count,
                                   self.vaf())

    # def corrected_depth(self):
    #     return int(self.depth) - self.bases['N'].count


def main():
    from pathlib import Path
    input_file = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    sample_file = output_path / f"{'.'.join(input_file.name.split('.')[0:1])}.csv"
    sample_id = '.'.join(sample_file.name.split('.')[0:1])
    sample_file = open(sample_file, 'w')

    for readcount in iter_readcounts(input_file):
        sample_position = (sample_id, # 1. library_id
                                readcount.position, # 2. position_id
                                readcount.allele_count(), # 3. num_alleles_observed
                                readcount.bases['A'].count, # 4. A_count
                                readcount.bases['C'].count, # 5. C_count
                                readcount.bases['G'].count, # 6. G_count
                                readcount.bases['T'].count, # 7. T_count
                                readcount.depth, # 8. allele_depth_count - does not include deletions
                                readcount.total_plus_count, # 9. forward_strand_count
                                readcount.plus_ratio(), # 10. forward_strand_ratio - FSR in paper
                                readcount.total_minus_count, # 11. reverse_strand_count
                                readcount.minus_ratio(), # 12. reverse_strand_ratio
                                readcount.alt_base().base, # 13. alternative_allele
                                readcount.alt_base().count, # 14. alternative_allele_count
                                readcount.alt_base().num_plus_strand, # 15. alternative_allele_forward_strand_count
                                readcount.alt_plus_vaf(), # 16. alternative_allele_forward_strand_ratio
                                readcount.alt_base().num_minus_strand, # 17. alternative_allele_reverse_strand_count
                                readcount.alt_minus_vaf(), # 18. alternative_allele_reverse_strand_ratio
                                readcount.vaf(), # 19. alternative_allele_frequency - AAF in paper
                                readcount.p_value()) # 20. strand_bias_likelihood - S in paper

        line = ','.join(map(str, sample_position))
        sample_file.write(line + '\n')
        
    sample_file.close()

from pprint import pprint

def test():
    test_path=""
    t = open(f'test_25511.base').read().strip()
    r = parse_line(t)

    from pathlib import Path

    input_file = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    sample_file = output_path / f"{input_file.name.split('.')[0]}.csv"
    pprint(sample_file)
    sample_id = '.'.join(sample_file.name.split('.')[0:1])
    pprint(sample_id)

    pprint(r)
    pprint(r.position)
    pprint(r.reference_base)
    pprint(r.alt_base())
    pprint(r.vaf())
    pprint(r.bases)
    pprint(r.plus_ratio())
    pprint(r.minus_ratio())
    pprint(r.alt_plus_vaf())
    pprint(r.alt_minus_vaf())
    pprint(r.depth)
    pprint(r.total_plus_count)
    pprint(r.total_minus_count)
    pprint(f"p-value is: {r.p_value()}")
    pprint('\n')

if __name__ == '__main__':
    #test()
    main()


