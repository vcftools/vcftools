<body>
<div id="wrap">
  <div id="top">
    <h2> <a href="index.html"><strong>VCFtools</strong></a></h2>
    <div id="menu">
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><a href="http://sourceforge.net/projects/vcftools/">Sourceforge page</a></li>	
        <li><a href="examples.html">Examples &amp; Documentation</a></li>
        <li><a href="downloads.html">Downloads</a></li>
      </ul>
    </div>
  </div>
  <div id="content">
      <h2>The C++ executable module examples</h2>
		<p>This page provides usage examples for the executable module. Extended documentation <b>for all of the options</b> can be found on the <a href="man_latest.html">manual page</a>. 
		</p>
		<ul>
			<li><a href="#run">Running the program</a></li>
			<li><a href="#file">Getting basic file statistics</a></li>
			<li><a href="#filter">Applying a filter</a></li>	
			<li><a href="#newvcf">Writing to a new VCF file</a></li>
			<li><a href="#streamvcf">Writing out to screen</a></li>
			<li><a href="#convertbcf">Converting a VCF file to BCF</a></li>
			<li><a href="#comparevcf">Comparing two VCF files</a></li>
			<li><a href="#freq">Getting allele frequency</a></li>
			<li><a href="#depth">Getting sequencing depth information</a></li>
			<li><a href="#ld">Getting linkage disequilibrium statistics</a></li>
			<li><a href="#fst">Getting Fst population statistics</a></li>
			<li><a href="#plink">Converting VCF files to PLINK format</a></li>
		</ul>
		
		<h2><a name="run" class="Q">Run the program</a></h2>
		<p>By default the executable can be found in the <b>bin/</b> subdirectory. To run the program, type:</p>
		<p class="codebox">./vcftools</p>
		<p>The program will return information regarding the version number.</p>
		    
		<h2><a name="file" class="Q">Get basic file statistics</a></h2>
		<p>The executable can be run with only an input VCF file without any other options, and will return basic information regarding the contents of the file.
		To specify an input file you must use the one of the input options ( <b>--vcf</b>, <b>--gzvcf</b>, or <b>--bcf </b> ) depending on the type of file.
		For example, for a VCF file called <b>input_data.vcf</b> the following command could be run:</p>
		<p class="codebox">./vcftools --vcf input_data.vcf</p>
		<p>It will return information about the file such as the <b>number of variants</b> and the <b>number of individuals</b> in the file.</p>
		<p>Beginning with vcftools v0.1.12, the program can also take input in from standard input (<b>stdin</b>). To do this, use any of the normal file type input options followed by the dash <b>-</b> character.</p>
		<p class="codebox">zcat input_data.vcf.gz | ./vcftools --vcf -</p>

		<h2><a name="filter" class="Q">Applying a filter</a></h2>
		<p>You can use VCFtools to filter out variants or individuals based on the values within the file.
		For example, to filter the sites within a file based upon their location in genome, use the options
		<b>--chr</b>, <b>--from-bp</b>, and <b>--to-bp</b> to specify the region.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --chr 1 --from-bp 1000000 --to-bp 2000000</p>
		<p>After running this line, the program will return the amount of sites in the file that are included in the chromosomal region chr1:1000000-2000000.
		This option can be modified to work with any desired region.</p>

		<h2><a name="newvcf" class="Q">Writing to a new VCF file</a></h2>
		<p>VCFtools can perform analyses on the variants that pass through the filters or simply <b>write those variants out to a new file</b>.
		This function is helpful for creating subsets of VCF files or just removing unwanted variants from VCF files. To write out the variants that pass through filters use the <b>--recode</b> option.
		In addition, use <b>--recode-INFO-all</b> to include all data from the INFO fields in the output. By default INFO fields are not written because many filters will alter the variants in a file, rendering
		the INFO values incorrect.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode --recode-INFO-all</p>
		<p>In this example, VCFtools will create a new VCF file containing only variants within the specified chromosomal region while keeping all INFO fields included in the original file.</p>
		<p>Any files written out by VCFtools will be in the current working directory and have the prefix <i>./out.SUFFIX</i> by default.
		To change the path, specify the new path using the option <b>--out</b> followed by the desired path. The program will add a suffix to that path based on the chosen output function.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode --out subset</p>

		<h2><a name="streamvcf" class="Q">Writing out to screen</a></h2>
		<p>Beginning with VCFtools v0.1.12, the program can also write out to screen instead of having the program write to a specified path. Using the options
		<b>--stdout</b> or <b>-c</b> will redirect all output to standard out. The output can then be piped into other programs or written out to a specified file name.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode --stdout | more</p>
		<p>The above example will output the resulting file to screen one line at a time for <b>quick inspection</b> of the results.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode -c > /home/usr/data/subset.vcf</p>
		<p>The above example will redirect the output and <b>write</b> it to the specified file name.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --chr 1 --from-bp 1000000 --to-bp 2000000 --recode -c | gzip -c > /home/usr/data/subset.vcf.gz</p>
		<p>The above example will redirect the output into <b>gzip</b> (assuming it is installed) for <b>compression</b>, and then gzip will write the file to the specified destination.</p>

		<h2><a name="convertbcf" class="Q">Converting a VCF file to BCF</a></h2>
		<p>Beginning with VCFftools v0.1.11, the program has the ability to read and write <b>BCF files</b>. This means that the program can also convert files between the two formats.
		This is accomplished in a similar way as the above example, instead using the <b>--recode-bcf</b> option. All output BCF files are automatically compressed using BGZF.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --recode-bcf --recode-INFO-all --out converted_output</p>
		
		<h2><a name="comparevcf" class="Q">Comparing two files</a></h2>
		<p>Using VCFtools, two VCF files can be compared to determine which sites and individuals are shared between them. The first file is declared using the input file options just like any other output function. The second file must be
		specified using <b>--diff</b>, <b>--gzdiff</b>, or <b>--diff-bcf</b>. There are also advanced options to determine additional discordance between the two files.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --diff other_data.vcf --out compare</p>
		
		<h2><a name="freq" class="Q">Getting allele frequency</a></h2>
		<p>To determine the frequency of each allele over all individuals in a VCF file, the <b>--freq</b> argument is used.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --freq --out output</p>
		<p>The output file will be written to <i>output.frq</i>.</p>
		
		<h2><a name="depth" class="Q">Getting sequencing depth information</a></h2>
		<p>Another useful output function <b>summarizes sequencing depth</b> for each individual or for each site. Just like the allele frequency example above, this output function follows the same basic model.
		</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --depth -c > depth_summary.txt</p>
		<p>With VCFtools, you can use many combinations of filters and an output function. For example, to write out site-wise sequence depths only at sites that have no missing data,
		include the <b>--max-missing</b> argument.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --site-depth --max-missing 1.0 --out site_depth_summary</p>		
		
		<h2><a name="ld" class="Q">Getting linkage disequilibrium statistics</a></h2>
		<p>Linkage disequilibrium between sites can be determined as well. This is accomplished using the <b>--hap-r2</b>, <b>--geno-r2</b>, or <b>--geno-chisq</b> arguments.
		Since the program must do pairwise site comparisons, this analysis can be time consuming, so it is recommended to filter the sites first or use one of the other options (<b>--ld-window</b>, <b>--ld-window-bp</b> or <b>--min-r2</b>) to reduce the number of comparisons.
		In this example, the VCFtools will only compare sites within 50,000 base pairs of one another.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --hap-r2 --ld-window-bp 50000 --out ld_window_50000</p>
		
		<h2><a name="fst" class="Q">Getting Fst population statistics</a></h2>
		<p>VCFtools can also calculate <b>Fst statistics</b> between individuals of different populations. It is an estimate calculated in accordance to Weir and Cockerham’s 1984 paper. The user must supply text files that contain <b>lists of individuals</b> (one per line) that are members of each population. The function will work with multiple populations if multiple <b>--weir-fst-pop</b> arguments are used.
		The following example shows how to calculate a per-site Fst calculation with two populations. Other arguments can be used in conjunction with this function, such as <b>--fst-window-size</b> and <b>--fst-window-step</b>.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --weir-fst-pop population_1.txt --weir-fst-pop population_2.txt --out pop1_vs_pop2</p>
		
		<h2><a name="plink" class="Q">Converting VCF files to PLINK format</a></h2>
		<p>VCFtools can convert VCF files into formats convenient for use in other programs. One such example is the ability to convert into PLINK format.
		The following function will output the variants in <b>.ped</b> and <b>.map</b> files.</p>
		<p class="codebox">./vcftools --vcf input_data.vcf --plink --chr 1 --out output_in_plink</p>
		
    <div id="clear"></div>
  </div>
  <div id="footer">
    <p>Copyright 2014 &copy; VCFtools</p>
  </div>
</div>
</body>
</html>
