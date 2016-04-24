<html>
<head>
<title>Automated Variations of Shower Parameters</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='Variations.php'>
 
<h2>Automated Variations of Shower Parameters for Uncertainty Bands</h2> 
<p> 
While a number of different central "tunes" of the Pythia parameters 
are provided, it is often desired  to study how event properties change when 
some of the parameters (such as those describing the parton showers) are 
varied.   Pythia8 now has the ability to provide a series of weights 
to reflect the change in probability for a particular final state to occur 
when a subset of parton-shower parameters are varied.  Details on the 
implementation and interpretation of these weights can be found in 
[<a href="Bibliography.php" target="page">Mre16</a>]. 
Currently, the list of available automated variations 
(see <a href="#keywords">full list below</a>) includes: 
<ul> 
<li> The renormalization scale for QCD emissions in FSR; </li> 
<li> The renormalization scale for QCD emissions in ISR; </li> 
<li> The inclusion of non-singular terms in QCD emissions in FSR; </li> 
<li> The inclusion of non-singular terms in QCD emissions in ISR. </li> 
</ul> 
Similar variations would be possible for QED emissions, but these have not 
yet been implemented. 
</p> 
 
<p> 
Since the computation of the uncertainty variations takes additional 
CPU time (albeit much less than would be required for independent runs 
with the equivalent variations), the automated uncertainty variations are 
switched off by default. They can be enabled by using the switch: 
<br/><br/><strong>UncertaintyBands:doVariations</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
   
</p> 
<h3>Specifying the Variations</h3> 
<p>When <code>UncertaintyBands:doVariations</code> is switched on, the user 
can define an arbitrary number of (combinations of) uncertainty variations 
to perform. Each variation is defined by a string with the following 
generic format:</p> 
<p> 
<code>Label keyword1=value keyword2=value ...</code> 
</p> 
<p>where the user has complete freedom to specify the label, and each 
keyword must be selected from the <a href="#keywords">list of currently 
recognised keywords below</a>. 
</p> 
<p>For example, an 
uncertainty variation corresponding to simultaneously increasing both 
the ISR and FSR renormalisation scales by a factor of two would be 
defined as follows 
</p> 
<p> 
<code>myVariation1 fsr:muRfac=2.0 isr:muRfac=2.0</code> 
</p> 
<p>Staying within the context of this example, the user might also want to 
check what a variation of the two scales independently of each other would 
produce. This can be achieved within the same run by adding two further 
variations, as follows: 
</p> 
<p> 
  <code>myVariation2 fsr:muRfac=2.0</code></br> 
  <code>myVariation3 isr:muRfac=2.0</code> 
</p> 
<p> 
Different histograms can then be filled with each set of weights as 
desired (see <a href="#access">accessing the uncertainty weights</a> below). 
  Variations by smaller or larger factors can obviously also be added in the 
  same way, again within one and the same run. 
</p> 
<p>Once a list of variations defined as above has been decided on, 
the whole list should be passed to Pythia in the form of a single 
"vector of strings", defined as follows:</p> 
<p/><code>wvec&nbsp; </code><strong> UncertaintyBands:List &nbsp;</strong> 
 (<code>default = <strong>{alphaShi fsr:muRfac=0.5 isr:muRfac=0.5,                       alphaSlo fsr:muRfac=2.0 isr:muRfac=2.0,                       hardHi fsr:cNS=2.0 isr:cNS=2.0,                       hardLo fsr:cNS=-2.0 isr:cNS=-2.0}</strong></code>)<br/>
  Vector of uncertainty-variation strings defining which variations 
  will be calculated by Pythia, when 
  <code>UncertaintyBands:doVariations</code> is switched on. 
   
<p>For completeness, we note that a command-file specification 
equivalent to the above default variations could look as follows:</p> 
<pre> 
UncertaintyBands:List = { 
  alphaShi fsr:muRfac=0.5 isr:muRfac=0.5, 
  alphaSlo fsr:muRfac=2.0 isr:muRfac=2.0, 
  hardHi fsr:cNS=2.0 isr:cNS=2.0, 
  hardLo fsr:cNS=-2.0 isr:cNS=-2.0 
} 
</pre> 
<p>Note that each of the individual uncertainty-variation definitions 
(the elements of the vector) are separated by commas and that 
keywords separated only by spaces are interpreted as belonging to a 
single combined variation. Note also that the beginning and end of the 
vector is marked by curly braces. 
</p> 
 
<a name="access"></a> 
<h3>Accessing the Uncertainty Weights</h3> 
<p>During the event generation, uncertainty weights will be calculated 
for each variation defined above, via the method described in 
[<a href="Bibliography.php" target="page">Mre16</a>]. 
The resulting alternative weights for the event are accessible through the 
   <code>Pythia::info.weight(int iWeight=0)</code> 
   method. 
</p> 
<p>The baseline weight for each event (normally unity for an 
ordinary unweighted event sample) is not modified and 
corresponds to <code>iWeight = 0</code>. The uncertainty-variation 
weights are thus enumerated starting from <code>iWeight = 1</code> for 
the first variation up to <code>N</code> for the last variation, in 
the order they were specified in <code>UncertaintyBands:List</code>. 
</p> 
<p>The total number of variations that have been defined, 
<code>N</code>, can be queried using 
<code>Pythia::info.nWeights()</code>. 
</p> 
<h3>NLO Compensation Term for Renormalisation-Scale Variations</h3> 
<p> 
  Additionally, there is a run-time parameter: 
 
<br/><br/><strong>UncertaintyBands:muSoftCorr</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
This flags tells the shower to apply an O(&alpha;S<sup>2</sup>) 
compensation term to the renormalization-scale variations, which 
reduces their magnitude for soft emissions, as described in 
[<a href="Bibliography.php" target="page">Mre16</a>]. 
   
</p> 
 
<a name="keywords"></a> 
<h3>List of Recognised Keywords for Uncertainty Variations</h3> 
<p> 
  The following keywords adjust the renormalisation scales and 
non-singular terms for all FSR and ISR branchings, respectively: 
<ul> 
 <li><code>fsr:muRfac</code> : multiplicative factor applied to the 
  renormalization scale for FSR branchings.</li> 
  <li><code>isr:muRfac</code> : multiplicative factor applied to the 
  renormalization scale for ISR branchings.</li> 
  <li><code>fsr:cNS</code> : additive non-singular ("finite") 
   term in the FSR splitting functions.</li> 
  <li><code>isr:cNS</code> : additive non-singular ("finite") 
   term in the ISR splitting functions.</li> 
</ul> 
Note that the <code>muRfac</code> parameters are applied linearly to 
the renormalisation scale, hence &mu;<sup>2</sup> &rarr; 
(<code>muRfac</code>)<sup>2</sup>*&mu;<sup>2</sup>;. 
</p> 
<p>Optionally, a further level of detail can be accessed by specifying 
variations for specific types of branchings, with the global keywords 
above corresponding to setting the same value for all 
branchings. Using the <code>fsr:muRfac</code> parameter for 
illustration, the individual branching types that can be specified 
are: 
<ul> 
<li><code>fsr:G2GG:muRfac</code> : variation for g&rarr;gg branchings.</li> 
<li><code>fsr:Q2QG:muRfac</code> : variation for q&rarr;qg branchings.</li> 
<li><code>fsr:G2QQ:muRfac</code> : variation for g&rarr;qqbar branchings.</li> 
<li><code>fsr:X2XG:muRfac</code> : variation for gluon bremsstrahlung off 
other types of particles (such as coloured new-physics particles). </li> 
</ul> 
For the distinction between <code>Q2QG</code> and <code>X2XG</code>, 
the following switch can be used to control whether b and t quarks are 
considered to be <code>Q</code> or <code>X</code> particles (e.g. providing a 
simple way to control top-quark or bottom-quark radiation 
independently of the rest of the shower uncertainties): 
<p/><code>mode&nbsp; </code><strong> UncertaintyBands:nFlavQ &nbsp;</strong> 
 (<code>default = <strong>6</strong></code>; <code>minimum = 2</code>; <code>maximum = 6</code>)<br/>
Number of quark flavours controlled via 
<code>Q2QG</code> keywords, with higher ID codes controlled by 
<code>X2XG</code> keywords. I.e., a change to 5 would mean that 
top-quark variations would use <code>X2XG</code> keyword values instead 
of the corresponding <code>Q2QG</code> ones. 
</modeselect> 
</p> 
<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "off")
{
$data = "UncertaintyBands:doVariations = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "UncertaintyBands:muSoftCorr = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2016 Torbjorn Sjostrand --> 
