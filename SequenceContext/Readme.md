<H3>Detection and enrichment of motifs</H3>

We collected for clusters A, E, G, H, L and M all SSMs of the subtype that is the most characteristic. This is C>A for clusters A and H, C>G for cluster E, C>T for cluster G, T>G for cluster L and C>G for cluster M. In addition, we looked at T>G SSMs in cluster H to compare to cluster L. Next, we extracted from the reference genome (GRCh37/h19) the ten adjacent bases in 5’ and 3’ direction of the mutation using the Rsamtools package. We used the extracted sequence context as input to construct two sequence logos per cluster: one for the mutations that are recurrent within the cluster and one for those that are not. As a measure of information content we used the relative entropy [50, 51], which is defined for position <em>i</em> by:

RE<sub>i</sub>=∑_<sub>b ∈{A,C,G,T}</sub>)f(b<sub>i</sub>) log_2  (f(b_i))/(P(b))

RE<sub>i</sub>=∑<sub>b∈{A,C,G,T}</sub>f(b<sub>i</sub>) log<sub>2</sub> &frasf(b>sub>i</sub>);<sub>P(b)</sub>
\frac{f(b>sub>i</sub>)}{P(b)}

<p>
    <span>&Sigma;</span>
    f(b<sub>i</sub>) log<sub>2</sub>  (f(b_i))/(P(b))
</p>

Here, f(b_i) stands for the frequency of base b (A, C, G or T) in position i and P(b) stands for the prior probability of base b as determined by the frequency in the human genome (GRCh37/h19). The height of each base in the sequence plot is proportional to f(b_i ) 〖log〗_2  (f(b_i))/(P(b)). A positive value corresponds to an enrichment of the base with respect to the prior probability and a negative value to a depletion. The relative entropy (REi) is zero, if all four bases are observed with the same frequency as the prior in position i. We set 0.25 as a threshold for REi to define the enriched motif. Furthermore, we computed per cluster the percentages of all, non-recurrent and recurrent SSMs that were in the sequence context that was found to be enriched in the recurrent SSMs. To estimate the percentage of the respective motifs in the human genome, we first slid a window of the same size (k) as the motif across the genome with a shift equal to the length of the motif and counted all possible k-mers. Next, we added to this the counts retrieved in the same way for the reverse complement of the reference sequence (corresponding to the opposing strand), since we also combined the reverse complements for each of the SSM subtypes. From this we computed the percentage of the enriched motif with respect to all k-mers and to the k-mer with the base that is mutated in the enriched motif at the same position.

<br>
<br>
50.	Schneider TD, Stormo GD, Gold L, Ehrenfeucht A. Information content of binding sites on nucleotide sequences. J Mol Biol. 1986;188(3):415-31. doi: 10.1016/0022-2836(86)90165-8.<br>
51.	Kullback S, Leibler RA. On Information and Sufficiency. Ann Math Statist. 1951;22(1):79-86.
