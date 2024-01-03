<?php


$data = $_REQUEST;

include('loading.html');

?>

<script>
	var dataObj = <?= json_encode($data) ?>;

	var data = { 	var1: dataObj.sequence,
			var2: dataObj.blast,
			var3: dataObj.lower,
			var4: dataObj.upper,
			var5: dataObj.pam_seq,
			var6: dataObj.sgrna,
			var7: dataObj.slen,
			var8: dataObj.mis,
			var9: dataObj.tmis
	 };
	 //console.log(dataObj);
	 
	window.location.href = "run_crispr.php?var1=" + dataObj.sequence + "&var2=" + dataObj.blast + "&var3=" + dataObj.lower + "&var4=" + dataObj.upper + "&var5=" + dataObj.pam_seq + "&var6=" + dataObj.sgrna + "&var7=" + dataObj.slen + "&var8=" + dataObj.mis + "&var9=" + dataObj.tmis;

</script>