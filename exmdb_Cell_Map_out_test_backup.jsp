<%@ page language="java" import="java.sql.*,wp.base.*" contentType="text/html; charset=UTF-8"
    pageEncoding="UTF-8"%>
<!DOCTYPE html>
<html>
<head>
<meta charset="ISO-8859-1">



<title>CeRNA-Cell-Map</title>
<!-- Main style sheet -->

<link rel="stylesheet" href="assets/css/bootstrap.min.css">
<link rel="stylesheet" href="assets/css/animate.min.css">
<link rel="stylesheet" href="assets/fonts/flaticon.css">
<link rel="stylesheet" href="assets/css/boxicons.min.css">
<link rel="stylesheet" href="assets/css/owl.carousel.min.css">
<link rel="stylesheet" href="assets/css/owl.theme.default.min.css">
<link rel="stylesheet" href="assets/css/magnific-popup.css">
<link rel="stylesheet" href="assets/css/nice-select.min.css">
<link rel="stylesheet" href="assets/css/meanmenu.css">
<link rel="stylesheet" href="assets/css/style.css">
<link rel="stylesheet" href="assets/css/responsive.css">
<link rel="icon" type="image/png" href="assets/images/favicon.png">

<link rel="stylesheet" type="text/css" href="select/demo.css"/>
<link rel="stylesheet" type="text/css" href="select/style-adsila.css" />
<link rel="stylesheet" href="select/selectpage.css" type="text/css">
<link rel="stylesheet" href="select/selectpage_cell.css" type="text/css">
<!-- echarts -->
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5.3.1/dist/echarts.min.js"></script>

<style>
.imm_select{
	width:30%;
}

.graph_qi {
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    z-index: 99999;
    background: #00000073;
    display:none;
}
</style>

<% 
String location="umap";
if(request.getParameter("location")!=null&&!request.getParameter("location").equals("")){	
	location = request.getParameter("location");
}

String resolution="5";
if(request.getParameter("resolution")!=null&&!request.getParameter("resolution").equals("")){	
	resolution = request.getParameter("resolution");
}

String character="celltype";
if(request.getParameter("character")!=null&&!request.getParameter("character").equals("")){	
	character = request.getParameter("character");
}

String gseid="GSE83142";
if(request.getParameter("character")!=null&&!request.getParameter("character").equals("")){	
	character = request.getParameter("character");
}
%>

</head>

<body>

<div class="graph_qi" name="graph_qi" id="graph_qi">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>

<div class="row">

<div class="col-lg-3 col-md-3">
<select class="language-list-item imm_select" name="gseid" id="gseid" onchange="change()">
		<option value="GSE83142">GSE83142</option>
		<option value="GSE142213" selected="selected">GSE142213</option>
		<option value="GSE113660">GSE113660</option>
</select>
</div>
<div class="col-lg-3 col-md-3">
<select class="language-list-item imm_select" name="location" id="location" onchange="change()">
		<option value="umap" selected="selected">umap</option>
		<option value="tsne">tsne</option>
</select>
</div>
<div class="col-lg-3 col-md-3">
<select class="language-list-item imm_select" name="character" id="character" onchange="change()">
		<option value="celltype">Cell Type</option>
		<option value="source">Source</option>
		<option value="stage">Stage</option>
		<option value="sample" selected="selected">Sample</option>
</select>
</div>
<div class="col-lg-3 col-md-3">
<select class="language-list-item imm_select" name="resolution" id="resolution" onchange="change()">
		<option value="1">0.1</option>
		<option value="2">0.2</option>
		<option value="3">0.3</option>
		<option value="4">0.4</option>
		<option value="5" selected="selected">0.5</option>
		<option value="6">0.6</option>
		<option value="7">0.7</option>
		<option value="8">0.8</option>
		<option value="9">0.9</option>
</select>
</div>

</div>







							<iframe frameborder=0 width="100%" height="500px" 
							name="exmdb_cell" id="exmdb_cell"
							src="exmdb_Cell_Map_res.jsp"
							scrolling="no" onload="load()"></iframe>
							
							<iframe frameborder=0 width="100%" height="500px" 
							name="exmdb_cell_time" id="exmdb_cell_time"
							src="exmdb_Cell_Map_time_res.jsp"
							scrolling="no" onload="load()"></iframe>
							
							<iframe frameborder=0 width="100%" height="1000px" 
							name="exmdb_cell" id="exmdb_cell"
							src="exmdb_singlecell_table_in.jsp"
							scrolling="no" onload="load()"></iframe>

<!-- left --> 

<script src="assets/js/jquery.min.js"></script>
<script src="assets/js/bootstrap.bundle.min.js"></script>
<script src="assets/js/owl.carousel.min.js"></script>
<script src="assets/js/jquery.magnific-popup.min.js"></script>
<script src="assets/js/jquery.nice-select.min.js"></script>
<script src="assets/js/wow.min.js"></script>
<script src="assets/js/meanmenu.js"></script>
<script src="assets/js/jquery.ajaxchimp.min.js"></script>
<script src="assets/js/form-validator.min.js"></script>
<script src="assets/js/contact-form-script.js"></script>
<script src="assets/js/custom.js"></script>

<script type="text/javascript" src="select/demo2.js" ></script>
<script type="text/javascript" src="select/jquery.min.js" ></script>
<script type="text/javascript" src="select/selectpage.min.js" ></script>

<script>
function change()
{
	document.getElementById("exmdb_cell").src="exmdb_Cell_Map_res.jsp?location="+document.getElementById("location").value+"&gseid="+document.getElementById("gseid").value+"&character="+document.getElementById("character").value+"&resolution="+document.getElementById("resolution").value;
	$("#graph_qi").fadeIn(1);
}



function load(){
	$("#graph_qi").fadeOut(1);
}
</script>

</body>
</html>