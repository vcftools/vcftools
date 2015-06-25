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
      <h2>The Perl modules examples</h2>
		<p>This page provides usage examples for the Perl modules. Extended documentation <b>for all of the options</b> can be found in the <a href="perl_module.html">full documentation</a>. 

		<ul>
			<li> <a href="#install">Installation</a>
			<li> <a href="#annotate">Annotating</a>
			<li> <a href="#compare">Comparing</a>
			<li> <a href="#concat">Concatenating</a>
			<li> <a href="#convert">Converting</a>
			<li> <a href="#filter">Filtering</a>
			<li> <a href="#isec">Intersections, complements</a>
			<li> <a href="#merge">Merging</a>
			<li> <a href="#query">Querying</a>
			<li> <a href="#shuffle">Reordering columns</a>
			<li> <a href="#stats">Stats</a>
			<li> <a href="#subset">Stripping columns</a>
			<li> <a href="#one-liners">Useful shell one-liners</a>
			<li> <a href="#validator">Validating</a>
		</ul>

		<h2><a name="annotate" class="Q">Annotating</a></h2>
		<p></p>
		<p class="codebox">
		# Add custom annotations <br>
		cat in.vcf | vcf-annotate -a annotations.gz \ <br>
		&nbsp;&nbsp; -d key=INFO,ID=ANN,Number=1,Type=Integer,Description='My custom annotation' \ <br>
		&nbsp;&nbsp; -c CHROM,FROM,TO,INFO/ANN &gt; out.vcf <br>
		<br>
		# Apply SnpCluster filter <br>
		cat in.vcf | <a href="perl_module.html#vcf-annotate">vcf-annotate</a> --filter SnpCluster=3,10 &gt; out.vcf <br>
		</p>

		<h2><a name="compare" class="Q">Comparing</a></h2>
		<p></p>
		<p class="codebox">
		<a href="perl_module.html#vcf-compare">vcf-compare</a> A.vcf.gz B.vcf.gz C.vcf.gz <br>
		<a href="htslib.html#check">vcf check</a> A.vcf.gz B.vcf.gz <br>
		</p>

		<h2><a name="concat" class="Q">Concatenating</a></h2>
		<p></p>
		<p class="codebox"><a href="perl_module.html#vcf-concat">vcf-concat</a> A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz
		</p>

		<h2><a name="convert" class="Q">Converting</a></h2>
		<p></p>
		<p class="codebox">
		# Convert between VCF versions  <br>
		zcat file.vcf.gz | <a href="perl_module.html#vcf-convert">vcf-convert</a> -r reference.fa | bgzip -c &gt; out.vcf.gz  <br>
		<br>
		# Convert from VCF format to tab-delimited text file <br>
		zcat file.vcf.gz | <a href="perl_module.html#vcf-to-tab">vcf-to-tab</a> &gt; out.tab <br>
		</p>

		<h2><a name="filter" class="Q">Filtering</a></h2>
		<p></p>
		<p class="codebox">
		# Filter by QUAL and minimum depth <br>
		<a href="perl_module.html#vcf-annotate">vcf-annotate</a> --filter Qual=10/MinDP=20 
		</p>

		<h2><a name="isec" class="Q">Intersections, complements</a></h2>
		<p></p>
		<p class="codebox">
		# Include positions which appear in at least two files  <br>
		<a href="perl_module.html#vcf-isec">vcf-isec</a> -o -n +2 A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz <br>
		<br>
		# Exclude from A positions which appear in B and/or C   <br>
		<a href="perl_module.html#vcf-isec">vcf-isec</a> -c A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz   <br>
		<br>
		# Fast hstlib implementation
		<a href="htslib.html#isec">vcf isec</a> -n =2 A.vcf.gz B.vcf.gz <br>
		</p>

		<h2><a name="merge" class="Q">Merging</a></h2>
		<p></p>
		<p class="codebox">
		<a href="perl_module.html#vcf-merge">vcf-merge</a> A.vcf.gz B.vcf.gz | bgzip -c &gt; C.vcf.gz <br>
		<a href="htslib.html#merge">vcf merge</a> A.vcf.gz B.vcf.gz <br>
		</p>

		<h2><a name="isec" class="Q">Querying</a></h2>
		<p></p>
		<p class="codebox"><a href="perl_module.html#vcf-query">vcf-query</a> file.vcf.gz 1:10327-10330 -c NA0001
		</p>

		<h2><a name="shuffle" class="Q">Reordering columns</a></h2>
		<p></p>
		<p class="codebox"><a href="perl_module.html#vcf-shuffle-cols">vcf-shuffle-cols</a> -t template.vcf.gz file.vcf.gz > out.vcf
		</p>

		<h2><a name="stats" class="Q">Stats</a></h2>
		<p></p>
		<p class="codebox"><a href="perl_module.html#vcf-stats">vcf-stats</a> file.vcf.gz <br>
		<a href="htslib.html#check">vcf check</a> file.vcf.gz &gt; file.vchk && plot-vcfcheck file.vchk -p plot/<br>
		</p>

		<h2><a name="subset" class="Q">Stripping columns</a></h2>
		<p></p>
		<p class="codebox"><a href="perl_module.html#vcf-subset">vcf-subset</a> -c NA0001,NA0002 file.vcf.gz | bgzip -c &gt; out.vcf.gz
		</p>

		<h2><a name="one-liners" class="Q">Useful shell one-liners</a></h2>
		<p> This sections lists some usefull one line commands. Note that there are also dedicated convenience scripts
		<a href="perl_module.html#vcf-sort">vcf-sort</a> and <a href="perl_module.html#vcf-concat">vcf-concat</a> 
		which do the same but also perform some basic sanity checks.
		All examples in BASH.</p>
		<p class="codebox">
		# Replace VCF header. The file must be compressed by bgzip. <br>
		tabix -r header.txt in.vcf.gz > out.vcf.gz                  <br>
		<br>
		# Sort VCF file keeping the header. The head command is for performance.    <br>
		(zcat file.vcf.gz | head -100 | grep ^#;                                    <br>
 		zcat file.vcf.gz | grep -v ^# | sort -k1,1d -k2,2n;)  \                    <br>
 		| bgzip -c &gt; out.vcf.gz <br>
		<br>
		# Merge (that is, concatenate) two VCF files into one, keeping the header   <br>
		#   from first one only.                    <br>
		(zcat A.vcf.gz | head -100 | grep ^#; \     <br>
 		zcat A.vcf.gz | grep -v ^#;          \     <br>
 		zcat B.vcf.gz | grep -v ^#; )        \     <br>
 		| bgzip -c &gt; out.vcf.gz     
		</p>

		<h2><a name="validator" class="Q">VCF validation</a></h2>
		<p> Both <span class="cmd">vcftools</span> and <span class="cmd">Vcf.pm</span> can be used for validation. The
    	first validates VCFv4.0, the latter is able to validate the older versions as well.
		</p>
		<p class="codebox">
		perl -M<a href="perl_module.html#Vcf.pm">Vcf</a> -e validate example.vcf  <br>
		perl -I/path/to/the/module/ -M<a href="perl_module.html#Vcf.pm">Vcf</a> -e validate example.vcf   <br>
		<a href="perl_module.html#vcf-validator">vcf-validator</a> example.vcf
		</p>

    <div id="clear"></div>
  </div>
  <div id="footer">
    <p>Copyright 2014 &copy; VCFtools</p>
  </div>
</div>
</body>
</html>
