<?php

function validate_path($path)
{
    if ( preg_match("/^[A-Za-z0-9_]+$/",$path) && file_exists("src/$path.inc") ) return $path;
    return "index";
}

$titles = array(
        'index'         => 'VCFtools',
        'perl_module'   => 'VCFtools: Perl tools and API',
        'perl_examples' => 'VCFtools: Perl tools and API Examples',
        'htslib'        => 'VCFtools: htslib VCF commands',
        'documentation' => 'VCFtools Documentation',
        'downloads'     => 'VCFtools Downloads',
        'examples'      => 'VCFtools Examples',
        'license'       => 'VCFtools License',
        'specs'         => 'VCF Specification',
        'links'         => 'VCF Links',
        'man_latest'		=> 'VCF Manual',
        );

if (isset($argc)) { $_GET['pg']=$argv[1]; }
$path  = array_key_exists('pg',$_GET) ? validate_path($_GET['pg']) : 'index';
$title = array_key_exists($path,$titles) ? $titles[$path] : $titles['index'];


?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
	<script type="text/javascript">
        var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
        document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
    </script>
    <script type="text/javascript">
        try {
            var pageTracker = _gat._getTracker("UA-272183-4");
            pageTracker._trackPageview();
        } catch(err) {}
    </script>
    
	<link href='favicon.png' rel='shortcut icon' type='image/png'>
  	<link href='favicon.png' rel='icon' type='image/png'>

<title><?php echo $title; ?></title>
<meta http-equiv="content-type" content="text/html; charset=UTF-8" />
<link rel="stylesheet" type="text/css" href="style.css" media="screen" />
</head>

<?php 
    include("$path.inc");
?>
<?php
