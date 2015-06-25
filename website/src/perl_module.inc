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

	<h2>The Perl modules and scripts</h2>
	<p>VCFtools contains a Perl API (<a href="#Vcf.pm">Vcf.pm</a>) and a number
	of Perl scripts that can be used to perform common tasks with VCF files such as
	file validation, file merging, intersecting, complements, etc. The Perl tools
	support all versions of the VCF specification (3.2, 3.3, 4.0, 4.1 and 4.2), 
	nevertheless, the users are encouraged to use the latest versions VCFv4.1 or VCFv4.2.
	The VCFtools in general have been used mainly with diploid data, but the Perl tools
	aim to support polyploid data as well.	Run any of the Perl scripts with the <b>--help</b> switch to obtain more help.	</p>
	
	<p>
	Many of the <b>Perl scripts require that the VCF files are compressed by <span class="cmd">bgzip</span>
	and indexed by <span class="cmd">tabix</span></b> (both tools are part of the tabix package,
	available for <a href="https://sourceforge.net/projects/samtools/files/tabix/">download here</a>).
	The VCF files can be compressed and indexed using the following commands
	<p class="codebox">bgzip my_file.vcf<br>
	tabix -p vcf my_file.vcf.gz
	</p>

	<h2>The tools</h2>
	<ul>
		<li> <a href="#fill-aa">fill-aa</a>
		<li> <a href="#fill-an-ac">fill-an-ac</a>
		<li> <a href="#fill-fs">fill-fs</a>
		<li> <a href="#fill-ref-md5">fill-ref-md5</a>
		<li> <a href="#fill-rsIDs">fill-rsIDs</a>
		<li> <a href="#vcf-annotate">vcf-annotate</a>
		<li> <a href="#vcf-compare">vcf-compare</a>
		<li> <a href="#vcf-concat">vcf-concat</a>
		<li> <a href="#vcf-consensus">vcf-consensus</a>
		<li> <a href="#vcf-contrast">vcf-contrast</a>
		<li> <a href="#vcf-convert">vcf-convert</a>
		<li> <a href="#vcf-filter">vcf-filter</a>
		<li> <a href="#vcf-fix-newlines">vcf-fix-newlines</a>
		<li> <a href="#vcf-fix-ploidy">vcf-fix-ploidy</a>
		<li> <a href="#vcf-indel-stats">vcf-indel-stats</a>
		<li> <a href="#vcf-isec">vcf-isec</a>
		<li> <a href="#vcf-merge">vcf-merge</a>
		<li> <a href="#vcf-phased-join">vcf-phased-join</a>
		<li> <a href="#vcf-query">vcf-query</a>
		<li> <a href="#vcf-shuffle-cols">vcf-shuffle-cols</a>
		<li> <a href="#vcf-sort">vcf-sort</a>
		<li> <a href="#vcf-stats">vcf-stats</a>
		<li> <a href="#vcf-subset">vcf-subset</a>
		<li> <a href="#vcf-to-tab">vcf-to-tab</a>
		<li> <a href="#vcf-tstv">vcf-tstv</a>
		<li> <a href="#vcf-validator">vcf-validator</a>
		<li> <a href="#Vcf.pm">Vcf.pm</a>
	</ul><br>

<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>
<script type="text/javascript">
jQuery(document).ready(function() {
        jQuery(".usageText").hide();
        jQuery(".usageToggle").click(function() { jQuery(this).next(".usageText").slideToggle(100); });
    });
</script>

<h3><a name="fill-an-ac" class="Q">fill-an-ac</a></h3>
<p> Fill or recalculate AN and AC INFO fields.
</p>
<p class="codebox">zcat file.vcf.gz | fill-an-ac | bgzip -c &gt; out.vcf.gz</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/fill-an-ac --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="fill-fs" class="Q">fill-fs</a></h3>
<p> Annotates the VCF file with flanking sequence (INFO/FS tag) masking known variants with N's. Useful for designing primers.
</p>
<p class="codebox">fill-fs -r /path/to/refseq.fa | vcf-query '%CHROM\t%POS\t%INFO/FS\n' &gt; out.tab</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/fill-fs --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="fill-ref-md5" class="Q">fill-ref-md5</a></h3>
<p>Fill missing reference info and sequence MD5s into VCF header.
</p>
<p class="codebox">fill-ref-md5 -i "SP:Homo\ Sapiens" -r ref.fasta in.vcf.gz -d ref.dict out.vcf.gz</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/fill-ref-md5 --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="fill-rsIDs" class="Q">fill-rsIDs</a></h3>
<p>Fill missing rsIDs. This script has been discontinued, please use <a href="#vcf-annotate">vcf-annotate</a> instead.
</p>

<h3><a name="vcf-annotate" class="Q">vcf-annotate</a></h3>
<p> The script adds or removes filters and custom annotations to VCF files.
To add custom annotations to VCF files, create TAB delimited file with annotations
such as
<p class="codebox" style="white-space:pre;font-family:monospace;">#CHR     FROM   TO      ANNOTATION 
1        12345  22345   gene1 
1        67890  77890   gene2 
</p>
<p>
Compress the file (using <span class="cmd">bgzip annotations</span>), index (using <span class="cmd">tabix -s 1 -b 2 -e 3 annotations.gz</span>) and
run
</p>
<p class="codebox">cat in.vcf | vcf-annotate -a annotations.gz \ <br>
&nbsp;&nbsp; -d key=INFO,ID=ANN,Number=1,Type=Integer,Description='My custom annotation' \ <br>
&nbsp;&nbsp; -c CHROM,FROM,TO,INFO/ANN &gt; out.vcf </p>
<p style='padding-bottom:0px'>The script is also routinely used to apply filters. There are a number of predefined filters and 
custom filters can be easily added, see <span class="cmd">vcf-annotate -h</span> for examples. 
Some of the predefined filters take advantage of tags added by bcftools, the descriptions of the most frequently asked
ones follow:</p>
    <div style='margin-left:1em;margin-top:0.2em;'>
        <b>Strand Bias</b> 
        .. Tests if variant bases tend to come from one strand. Fisher's exact test
        for 2x2 contingency table where the row variable is being the reference allele
        or not and the column variable is strand. Two-tail P-value is used.
    </div>
    <div style='margin-left:1em;margin-top:0.2em;'>
        <b>End Distance Bias</b> 
        .. Tests if variant bases tend to occur
        at a fixed distance from the end of reads, which is usually an indication of
        misalignment. (T-test)
    </div>
    <div style='margin-left:1em;margin-top:0.2em;'>
        <b>Base Quality Bias</b> 
        .. Tests if variant bases tend to occur with a quality bias (T-test). This filter is by default 
        effectively disabled as it is set to 0.
    </div>
<p><b><i>Note:</i></b> A fast <a href="htslib.html#annotate">htslib</a> C version of this tool is now available (see <a href="htslib.html#annotate">bcftools annotate</a>).
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-annotate --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
<div class="usageBox"><span class="usageToggle">(Read even more)</span><pre class="usageText"><?php unset($out);exec("cat ../examples/filters.txt",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-compare" class="Q">vcf-compare</a></h3>
<p> Compares positions in two or more VCF files and outputs the numbers of
positions contained in one but not the other files; two but not the other files, etc, which
comes handy when generating Venn diagrams. The script also computes numbers
such as nonreference discordance rates (including multiallelic sites), compares
actual sequence (useful when comparing indels), etc.
<br>
</p>
<p class="codebox">vcf-compare -H A.vcf.gz B.vcf.gz C.vcf.gz</p>
<p><b><i>Note:</i></b> A fast <a href="htslib.html#stats">htslib</a> C version of this tool is now available (see <a href="htslib.html#stats">bcftools stats</a>).
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-compare --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-concat" class="Q">vcf-concat</a></h3>
<p>Concatenates VCF files (for example split by chromosome). Note that the input and output VCFs will
have the same number of columns, the script does not merge VCFs by position (see also
<a href="#vcf-merge">vcf-merge</a>).
<p>In the basic mode it does not do anything fancy except for a sanity check that all
files have the same columns.  When run with the -s option, it will perform a partial
merge sort, looking at limited number of open files simultaneously.
<p class="codebox">vcf-concat A.vcf.gz B.vcf.gz C.vcf.gz | gzip -c &gt; out.vcf.gz</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-concat --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-consensus" class="Q">vcf-consensus</a></h3>
<p>Apply VCF variants to a fasta file to create consensus sequence.
</p>
<p class="codebox">cat ref.fa | vcf-consensus file.vcf.gz &gt; out.fa</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-consensus --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-convert" class="Q">vcf-convert</a></h3>
<p>Convert between VCF versions, currently from VCFv3.3 to VCFv4.0.
</p>
<p class="codebox">zcat file.vcf.gz | vcf-convert -r reference.fa &gt; out.vcf</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-convert --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-contrast" class="Q">vcf-contrast</a></h3>
<p> A tool for finding differences between groups of samples, useful in trio analysises, cancer genomes etc. 
<p>
In the example below variants with average mapping quality of 30 (<span class="cmd">-f MinMQ=30</span>) 
and minimum depth of 10 (<span class="cmd">-d 10</span>) are considered.
Only novel alleles are reported (<span class="cmd">-n</span>). Then <span class="cmd">vcf-query</span> 
is used to extract the INFO/NOVEL* annotations into a table. Finally the sites are 
sorted by confidence of the site being different in the child (<span class="cmd">-k5,5nr</span>).
</p>
<p class="codebox">
vcf-annotate -f MinMQ=30 file.vcf | vcf-contrast -n +Child -Mother,Father -d 10 -f | vcf-query -f '%CHROM %POS\t%INFO/NOVELTY\t%INFO/NOVELAL\t%INFO/NOVELGT[\t%SAMPLE %GTR %PL]\n' | sort -k3,3nr | head</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-contrast --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-filter" class="Q">vcf-filter</a></h3>
<p>Please take a look at <a href="#vcf-annotate">vcf-annotate</a> and <a href="htslib.html#view">bcftools view</a> which does what you are looking for. Apologies for the non-intuitive naming.
<br><b><i>Note:</i></b> A fast <a href="htslib.html">HTSlib</a> C version of a filtering tool is now available (see <a href="htslib.html#filter">bcftools filter</a> and
<a href=htslib.html#view>bcftools view</a>).
</p>

<h3><a name="vcf-fix-newlines" class="Q">vcf-fix-newlines</a></h3>
<p> Fixes diploid vs haploid genotypes on sex chromosomes, including the pseudoautosomal regions.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-fix-newlines --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-fix-ploidy" class="Q">vcf-fix-ploidy</a></h3>
<p> Fixes diploid vs haploid genotypes on sex chromosomes, including the pseudoautosomal regions.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-fix-ploidy --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-indel-stats" class="Q">vcf-indel-stats</a></h3>
<p> Calculate in-frame ratio.
</p>
<p><b><i>Note:</i></b> A fast <a href="htslib.html#stats">htslib</a> C version of this tool is now available (see <a href="htslib.html#stats">bcftools stats</a>).
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-indel-stats --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-isec" class="Q">vcf-isec</a></h3>
<p> Creates intersections and complements of two or more VCF files. Given multiple
    VCF files, it can output the list of positions which are shared by at least N files,
    at most N files, exactly N files, etc. The first example below outputs positions shared by
    at least two files and the second outputs positions present in the files A but absent
    from files B and C.
</p>
<p class="codebox">
vcf-isec -n +2 A.vcf.gz B.vcf.gz | bgzip -c &gt; out.vcf.gz <br>
vcf-isec -c A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz
</p>
<p><b><i>Note:</i></b> A fast <a href="htslib.html">htslib</a> C version of this tool is now available (see <a href="htslib.html#isec">bcftools isec</a>).
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-isec --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-merge" class="Q">vcf-merge</a></h3>
<p> Merges two or more VCF files into one so that, for example, if two source files
    had one column each, on output will be printed a file with two columns.
    See also <a href="vcf-concat">vcf-concat</a> for concatenating VCFs split by chromosome.
</p>
<p class="codebox">vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c &gt; out.vcf.gz</p>
<p>Note that this script is not intended for concatenating VCF files. For this, use
<a href="perl_module.html#vcf-concat">vcf-concat</a> instead.
<br><b><i>Note:</i></b> A fast <a href="htslib.html">htslib</a> C version of this tool is now available (see <a href="htslib.html#merge">bcftools merge</a>).
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-merge --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-phased-join" class="Q">vcf-phased-join</a></h3>
<p> Concatenates multiple overlapping VCFs preserving phasing.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-phased-join --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-query" class="Q">vcf-query</a></h3>
<p> Powerful tool for converting VCF files into format defined by the user.
Supports retrieval of subsets of positions, columns and fields.
</p>
<p class="codebox">
vcf-query file.vcf.gz 1:10327-10330<br>
vcf-query file.vcf -f '%CHROM:%POS %REF %ALT [ %DP]\n' 
</p>
<p>
<b><i>Note:</i></b> A fast <a href="htslib.html">htslib</a> C version of this tool is now available (see <a href="htslib.html#query">bcftools query</a>).
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-query --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="shuffle" class="Q">vcf-shuffle-cols</a></h3>
<p>Reorder columns</p>
<p class="codebox">vcf-shuffle-cols -t template.vcf.gz file.vcf.gz > out.vcf
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-shuffle-cols --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-sort" class="Q">vcf-sort</a></h3>
<p>Sort a VCF file.
</p>
<p class="codebox">vcf-sort file.vcf.gz</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-sort --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-stats" class="Q">vcf-stats</a></h3>
<p>Outputs some basic statistics: the number of SNPs, indels, etc.
</p>
<p class="codebox">vcf-stats file.vcf.gz</p>
<p><b><i>Note:</i></b> A fast <a href="htslib.html#stats">htslib</a> C version of this tool is now available (see <a href="htslib.html#stats">bcftools stats</a>).
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-stats --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-subset" class="Q">vcf-subset</a></h3>
<p>Remove some columns from the VCF file.
</p>
<p class="codebox">vcf-subset -c NA0001,NA0002 file.vcf.gz | bgzip -c &gt; out.vcf.gz</p>
<p>
<b><i>Note:</i></b> A fast <a href="htslib.html">HTSlib</a> C version of this tool is now available (see <a href="htslib.html#view">bcftools view</a>).
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-subset --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-tstv" class="Q">vcf-tstv</a></h3>
<p>A lightweight script for quick calculation of <a href="http://en.wikipedia.org/wiki/Transition_%28genetics%29">Ts</A>/<A href="http://en.wikipedia.org/wiki/Transversion">Tv</A> ratio.
</p>
<p class="codebox">cat file.vcf | vcf-tstv</p>
<p><b><i>Note:</i></b> A fast <a href="htslib.html#stats">htslib</a> C version of this tool is now available (see <a href="htslib.html#stats">bcftools stats</a>).
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-tstv --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-to-tab" class="Q">vcf-to-tab</a></h3>
<p>A simple script which converts the VCF file into a tab-delimited text file listing the
    actual variants instead of ALT indexes.
</p>
<p class="codebox">zcat file.vcf.gz | vcf-to-tab &gt; out.tab</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-to-tab --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="vcf-validator" class="Q">vcf-validator</a></h3>
<p>
</p>
<p class="codebox">vcf-validator file.vcf.gz</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("../perl/vcf-validator --help 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>

<h3><a name="Vcf.pm" class="Q">Vcf.pm</a></h3>
<p>For examples how to use the Perl API, it is best to look at some of the simpler scripts, 
for example <span class="cmd">vcf-to-tab</span>. The detailed documentation can be obtained by running</p>
<p class="codebox">perldoc Vcf.pm</p>

</div>
<div id="clear"></div>
</div>
<div id="footer">
    <p>Copyright 2014 &copy; VCFtools</p>
  </div>
</div>
</body></html>
