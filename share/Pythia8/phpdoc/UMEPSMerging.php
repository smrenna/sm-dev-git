 
<html>
<head>
<title>UMEPS Merging</title>
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

<form method='post' action='UMEPSMerging.php'>
 
<h2>Unitarised Matrix Element + Parton Shower Merging</h2> 
 
Pythia offers the possibility to use the unitarised matrix element + parton 
shower merging scheme, as presented in [<a href="Bibliography.php" target="page">Lon12</a>]. Unitarised ME+PS 
merging (UMEPS) allows for a consistent inclusion of tree-level multi-parton 
matrix elements into Pythia, and prevents potential changes in the inclusive 
production cross section. This makes it theoretically more appealing than 
CKKW-L merging. As in CKKW-L, UMEPS merging requires the user to supply Les 
Houches Event File input. 
 
<p/> 
UMEPS is different from other tree-level merging schemes in that it contains 
events with negative weights. These are generated by constructing 
parts of no-emission probabilities by reweighted higher-multiplicity 
samples [<a href="Bibliography.php" target="page">Lon12</a>]. The main philosophy of UMEPS is "subtract what you 
add", meaning that in order to ensure the stability of the inclusive cross 
section, one has to counter the inclusion of additional tree-level matrix 
elements by "subtraction terms". 
 
<p/> The scheme closely reflects how unitarity 
is achieved in a non-merged shower, and indeed explicitly enforces the 
cancellations that are implicitly happening in a non-merged shower. This makes 
very low merging scale values possible. 
 
<p/> 
The usage of UMEPS is illustrated in the sample main program 
<code>main86.cc</code>, together with the input file 
<code>main86.cmnd</code>. 
 
<p/> 
Unitarised merging is heavily indebted to CKKW-L merging, and shares many 
settings with CKKW-L. In particular, 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; The hard process 
(<code>Merging:Process</code>)needs to be defined 
exactly as in CKKW-L (see <strong>Defining the hard process</strong> in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L documentation</a>). 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; The merging scale value 
(<code>Merging:TMS</code>) has to be set. 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; The maximal number of additional partons 
<code>Merging:nJetMax</code> has to be set. 
 
<p/> 
UMEPS further shares the switches listed under the sections "<strong>Matrix 
element merging and HepMC output for RIVET</strong>" and "<strong>Further 
variables</strong>" in  the <?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L 
documentation</a> with CKKW-L merging. Also, all 
<code>MergingHooks</code> routines that allow for user interference in 
CKKW-L merging are also usable for UMEPS -- with the exception of a 
user-defined merging scale. 
Currently, UMEPS is only implemented for a merging scale defined 
by the minimal Pythia evolution pT value between sets of radiator, emitted 
and recoiler partons. This is no fundamental limitation of the method, and 
will possibly be lifted in the future. Since this merging scale definition is 
not completely obvious, UMEPS also shares the 
<code>Merging:enforceCutOnLHE</code> switch with CKKW-L. In this way, it 
is possible to use LHE files that are regularised only with weak cuts as 
input, while the merging machinery imposes the stronger merging scale cut 
automatically. This means that no merging scale implementation is required 
from the user side, but also means that it is the user's responsibility to 
ensure that the cuts used for generating input LHE files are always looser 
than the cut given by the merging scale value <code>Merging:TMS</code>. 
 
<br/><br/><hr/> 
<h3>UMEPS merging with main86.cc</h3> 
 
The UMEPS procedure is illustrated in the sample main program 
<code>main86.cc</code> (with the input card <code>main86.cmnd</code>). This 
program produces HepMC events [<a href="Bibliography.php" target="page">Dob01</a>], that can be histogrammed (e.g. 
using RIVET [<a href="Bibliography.php" target="page">Buc10</a>]), or used as input for a detector simulation. If 
the user is not familiar with HepMC analysis tools, it is possible to instead 
use Pythia's histogramming routines. For this, remove the lines referring to 
HepMC, and histogram events as illustrated (for CKKW-L) for the histogram 
<i>histPTFirstSum</i> in <code>main84.cc</code>, i.e. using 
<i>weight*normhepmc</i> as weight. 
 
<p/> 
In principle, no changes to <code>main86.cc</code> are necessary. Instead, all 
settings can be transferred to <code>main86.cc</code> through an input file. 
The input LHE files are also part of the (command line) input of 
<code>main86.cc</code>. Note  that the sample program assumes that LHE file 
names are of the form <i>name_tree_#nAdditionalJets.lhe</i>. If you want to 
e.g. use the LHE files that are shipped with the Pythia distribution, you can 
execute <code>main86.exe</code> with the command 
</p> 
<code>./main86.exe ./main86.cmnd ./w_production ./myhepmc.hepmc</code> 
</p> 
 
Since <code>main86.cc</code> is currently the "front-end" for UMEPS merging, 
we will briefly discuss this sample program in the following. 
 
<h4>Inputs</h4> 
 
In its current form, <code>main86.cc</code> uses separate tree-level LHE files 
for different numbers of additional partons as input. If e.g. UMEPS merging 
for W-boson + up to two additional partons is to be performed, three LHE files 
(for W+zero, W+one, W+two partons) are required. The configurations in the 
input files should be regularised with inclusive (i.e. weak) cuts. The actual 
"merging scale cut" will be handled internally. If e.g. 
<code>Merging:TMS = 15</code> is the desired merging scale value, 
it is acceptable 
to regularise the matrix element calculation for Higgs+jets events at the LHC 
with the loose cuts <i>pT<sub>jet</sub> = 5 GeV</i>, 
<i> &Delta;R<sub>jetA jetB</sub> = 0.01</i> and 
<i> Q<sub>jetA jetB</sub> = 5 GeV</i>. 
 
<p/> 
All input settings are handed to <code>main86.cc</code> in the form of an 
input file. This input file has to contain 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; The number of desired events 
(<code>Main:numberOfEvents</code>) 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; The hard process 
(<code>Merging:Process</code>) 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; The merging scale value 
(<code>Merging:TMS</code>) 
<p/> 
 &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; The maximal number of additional partons 
(<code>Merging:nJetMax</code>). 
 
<p/> 
Other settings are of course allowed. However, please refrain from adding 
switches that are used to invoke other merging schemes (e.g. 
<code>Merging:doKTMerging</code>) into the input file, since this can 
cause problems. 
 
 
<h4>Program flow</h4> 
 
The sample program starts by estimating the cross section for samples with 
different jet multiplicities. For this, the switch 
<code>Merging:doXSectionEstimate</code> is invoked together with the merging 
scale definition of <code>Merging:doUMEPSTree</code>, which corresponds to the 
minimal Pythia evolution pT value. We will come back to the latter switch 
below. All showering, multiparton interactions and hadronisation is, for speed 
reasons, switched off when estimating the cross section, since the hard cross 
section estimate would not be influenced by the event evolution anyway. 
 
<p/> 
After the hard cross sections are known (including the application of the 
merging scale cut), the first part of the UMEPS events is generated by 
using the following switch. 
 
<br/><br/><strong>Merging:doUMEPSTree</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Reweight events according to the UMEPS prescription for tree-level 
configurations. 
   
 
<p/> 
The weight generated by the UMEPS procedure can be accessed by using the 
function <strong> double Info::mergingWeight()</strong>. 
When printing (or histogramming) merged events, this weight, multiplied 
with the estimated cross section for the current sample, should be 
used as event weight (or to fill histogram bins). 
 
<p/> 
After this first part is complete, the outcome is an addition of reweighted 
tree-level samples. To restore the inclusive cross section (i.e. that the 
cross section after merging corresponds to the cross section of the hard 
process, without any additional jets), it is necessary to subtract samples. 
Parton shower unitarity leads to the conclusion that "resolved" and 
"unresolved" corrections always cancel between states that contain an 
additional resolved jet, and states in which we "integrate over" the phase 
space of the additional jet. 
<code>main86.cc</code> makes this cancellation explicit by producing 
(correctly weighted) counter events by switching on 
 
<br/><br/><strong>Merging:doUMEPSSubt</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Reweight events according to the UMEPS prescription of reweighted, 
integrated configurations. Please note that, in order for this to work 
smoothly, the switch <code>Merging:doUMEPSTree</code> has to be turned off. 
   
 
<p/> The integration is achieved internally, and the number of desired 
integrations (which is always one for UMEPS counter events) is set by 
 
<br/><br/><table><tr><td><strong>Merging:nRecluster  </td><td></td><td> <input type="text" name="3" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Number of hard partons to integrate out in the UMEPS procedure. 
   
 
<p/> Again, the weight generated by the UMEPS procedure can be accessed by 
using the function <strong> double Info::mergingWeight()</strong>. This 
weight, multiplied with the cross section of the current sample, and 
multiplied by <i>-1</i>, should then be used as event weight (or to fill 
histogram bins). 
 
<p/> 
Before returning, <code>main86.cc</code> prints the merged cross section 
after UMEPS merging. 
 
 
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
$data = "Merging:doUMEPSTree = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Merging:doUMEPSSubt = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0")
{
$data = "Merging:nRecluster = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2017 Torbjorn Sjostrand --> 
