<html>
<head>
<title>Variations of Shower Parameters</title>
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
 
<h2>Variations of Shower Parameters</h2> 
<p> 
While a number of different central Tunes are provided, it is often desired 
to study how event properties change when some of the parameters 
(such as those describing the parton showers) are 
varied.   Pythia8 now has the ability to provide a series of weights 
to reflect the change in probability for a particular final state to occur when 
a subset of parton shower parameters are varied.   Currently, the list of 
available variations are for: 
<ol> 
<li> the renormalization scale for QCD emissions in FSR; </li> 
<li> the renormalization scale for QCD emissions in ISR; </li> 
<li> the inclusion of finite pieces in QCD emissions in FSR; </li> 
<li> the inclusion of finite pieces in QCD emissions in ISR. </li> 
</ol> 
Similar variations are possible for QED emissions, but these are not 
yet implemented. 
</p> 
 
The variations are steered by a settings string with the format: 
 
<p>"user-provided-name-for-this-variation  variation-parameter value ... variation-parameter value"</p> 
 
Examples are: 
 
<ul> 
<li> "UncertaintyVariation A fsr:muR 0.5 isr:muR 0.5" </li> 
<li> "UncertaintyVariation C fsr:cNS 3.0" </li> 
</ul> 
 
<p/> 
The weights for each variation are available on an event-by-event basis 
using: 
<i> pythia.info.weight(int i); </i> 
where <i> i </i> represents the <i>i</i>'th variation. 
<i> weight(0) </i> is the baseline (unvaried) prediction with unit weight. 
 
<p/> 
 
Additionally, there is a run-time parameter: 
 
<br/><br/><strong>Uncertainties:muSoftCorr</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
This flags tells the shower to apply a soft correction to 
gluon emission. 
   
 
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
$data = "Uncertainties:muSoftCorr = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2016 Torbjorn Sjostrand --> 
