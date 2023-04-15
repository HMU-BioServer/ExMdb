<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
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

<link rel="stylesheet" href="bootstrap_select/bootstrap-select.min.css" type="text/css">

<title>ExMdb : Blast Analyze Tool</title>



<style>
.imm_select{
	width:100%;
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
.menu__item{
	margin: .5em 0;
}
</style>

</head>

<% 
String searchname = "HOTAIR";  
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
	searchname = request.getParameter("searchname");
}

String organ = "lung";  
if(request.getParameter("organ")!=null&&!request.getParameter("organ").equals("")){
	organ = request.getParameter("organ");
}


// String sql = "select * from exmdb_net where node1 in (\"A1BG-AS1\",\"MALAT1\",\"XIST\",\"miR-3605-5p\") or node2 in (\"A1BG-AS1\",\"MALAT1\",\"XIST\",\"miR-3605-5p\")";
// String data = dbhello.exmdb_netmaker(sql);
// System.out.println(data);
%>


<body>
<div class="preloader">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>

<div class="graph_qi" name="graph_qi" id="graph_qi">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>


<div class="navbar-area" data-step="1" data-intro="This is a tooltip!">
  <div class="mobile-nav"><a href="index.html" class="logo"><img src="assets/images/logos/logo-1.png" alt="Logo"></a></div>
  <div class="main-nav">
    <div class="container">
      <nav class="navbar navbar-expand-md navbar-light "><a class="navbar-brand" href="index.html"><img src="assets/images/logos/logo-1.png" alt="Logo"></a>
        <div class="collapse navbar-collapse mean-menu" id="navbarSupportedContent">
          <ul class="navbar-nav m-auto" style="z-index:100">
            <li class="nav-item"><a href="index.html" class="nav-link">Home</a></li>
            <li class="nav-item"><a href="exmdb_browse.jsp" target="_blank" class="nav-link">Browse</a></li>
            <li class="nav-item"><a href="#" class="nav-link">Search <i class='bx bx-caret-down'></i></a>
              <ul class="dropdown-menu">
                <li class="nav-item"><a href="exmdb_search_highput.jsp" class="nav-link" target="_blank">Search [ High throughput ]</a></li>
                <li class="nav-item"><a href="exmdb_search_exp.jsp" class="nav-link" target="_blank">Search [ Experimental validation ]</a></li>
                <li class="nav-item"><a href="exmdb_search_biomarker.jsp" class="nav-link" target="_blank">Search [ Cancer Biomarkers ]</a></li>
                <li class="nav-item"><a href="exmdb_search_scdataset.jsp" class="nav-link" target="_blank">Search [ Single-cell ]</a></li>
              </ul>
            </li>
            <li class="nav-item"><a href="#" class="nav-link active">Tools <i class='bx bx-caret-down'></i></a>
              <ul class="dropdown-menu">
                <li class="nav-item"><a href="exmdb_cell_out.jsp" target="_blank" class="nav-link">ExMdb-Cellloaction </a></li>
                <li class="nav-item"><a href="exmdb_imm_out.jsp" target="_blank" class="nav-link">ExMdb-Immunity </a></li><li class="nav-item"><a href="#" class="nav-link">ExMdb-Survival Tools<i class='bx bx-caret-down'></i></a><ul class="dropdown-menu"><li class="nav-item"><a href="exmdb_survival_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">mRNA-Survival Analyze</a></li><li class="nav-item"><a href="exmdb_survival_miRNA_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">miRNA-Survival Analyze</a></li> <li class="nav-item"><a href="exmdb_survival_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">lncRNA-Survival Analyze</a></li></ul></li>
                <li class="nav-item"><a href="exmdb_net_out.jsp" target="_blank" class="nav-link">ExMdb-Network </a></li><li class="nav-item"><a href="exmdb_blast_out.jsp" target="_blank" class="nav-link">ExMdb-BLAST </a></li>
                <li class="nav-item"><a href="exmdb_function_time.jsp" target="_blank" class="nav-link">ExMdb-Pseudotime pathway</a></li><li class="nav-item"><a href="exmdb_single_vis.jsp" target="_blank" class="nav-link">ExMdb-Single Cell Visualization</a></li>
				<li class="nav-item"><a href="#" class="nav-link">ExMdb-FunctionEnrichment <i class='bx bx-caret-down'></i></a><ul class="dropdown-menu"><li class="nav-item"><a href="exmdb_function_mRNA_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">mRNA-FunctionEnrichment</a></li><li class="nav-item"><a href="exmdb_function_out_miRNA.jsp" target="_blank" class="nav-link" style="text-transform: none;">miRNA-FunctionEnrichment</a></li> <li class="nav-item"><a href="exmdb_function_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">lncRNA-FunctionEnrichment</a></li></ul></li>
              </ul>
            </li>


            <li class="nav-item"><a href="exmdb_statistic.jsp" target="_blank" class="nav-link">Statistic </a></li><li class="nav-item"><a href="exmdb_help.html" target="_blank" class="nav-link">Help </a></li>
            
            <li class="nav-item"><a href="exmdb_download.html" target="_blank" class="nav-link">Download </a></li><li class="nav-item"><a href="exmdb_contact.html" target="_blank" class="nav-link">Contact Us</a></li>
          </ul>
          
        </div>
      </nav>
    </div>
  </div>
</div>

<div class="inner-banner">
  <div class="container">
    <div class="inner-title text-center">
      <h3>Blast Analyze Tool</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Blast Analyze Tool</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
<div class="terms-conditions-area" style="padding:20px 0px">
<div class="container wow fadeInUp animated" style="max-width: 1440px;">
		<div class="contact-header">
			<h1>Blast Analyze Tool</h1>
		</div>	
		<div class="contact-wrapper" style="padding: 1px 30px;">
				<form id="gene" name="gene" method="post" action="http://www.bio-server.cn/Server_blast/blast_res_new.jsp" target="exmdb_blast">
				
				<p class="sub_title_function" style="margin-top:20px"><img src="img/iconfont/blast.png" style="margin-top: -7px;margin-right: 7px;">ExMdb-BLAST</p>
	    		<p class="sub_text_function">The BLAST tool is convenient for users to query dataset by inputting custom sequences. Please input a sequence in the following textbox.</p>
				
				 <nav class="menu menu--adsila">
					<a class="menu__item" href="#">
						<span class="menu__item-name">ExMdb - BLAST</span>
						<span class="menu__item-label">( Gene sequence )</span>
					</a>

		        <hr>
				<div class="row">
		        <div class="col-md-10 col-sm-10 padd-0" style="padding:10px 30px;border-right: 1px solid #a0a0a0;">
		        	<textarea id="searchName" name="searchName" style="width:100%;height:100px;font-size:18px;padding: 10px;border: 2px solid #b3b3b3;" placeholder="Please input a sequence"></textarea>   
		        </div>
		        <div class="col-md-2 col-sm-2 padd-0">
		        	<button class="default-btn btn-bg-two border-radius-50" type="submit" style="width:100%;padding: 5px;" onclick="return gene_click()">BLAST</button><hr>
					<p style="font-size:9px">An example sequence of lncRNA "<a style="text-decoration: none; cursor: pointer;margin-left:0px;" onclick="document.getElementById('searchName').value='ACAGCAGGCAGCTGTTAACAGATAAGTTTAACTTGCATCTGCAGTATTGCATGTTAGGGATAAGTGCTTATTTTTAAGAGCTGTGGAGTTCTTAAATATCAACCATGGCACTTTCTCCTG'"><span style="color:#dcb200">MALAT1</span></a>".</p>
		        </div>
		        </div> 
				</nav>
				</form>

				
				<form id="gene1" name="gene1" method="post" action="http://www.bio-server.cn/Server_blast/blast_res_new1.jsp" target="exmdb_blast"> 
		        <!-- <form id="gene1" name="gene1" method="post" action="http://www.bio-bigdata.net/LncACTdb/blast_res_new.jsp" target="detail_res_pre_location"> -->
		        <nav class="menu menu--adsila">
					<a class="menu__item" href="#">
						<span class="menu__item-name">ExMdb - BLAST</span>
						<span class="menu__item-label">( Protein sequence )</span>
					</a>
				<hr>
				<div class="row">
		        <div class="col-md-10 col-sm-10 padd-0" style="padding:10px 30px;border-right: 1px solid #a0a0a0;">
		        	<textarea id="searchName1" name="searchName" style="width:100%;height:110px;font-size:18px;padding: 10px;border: 2px solid #b3b3b3;" placeholder="Please input a sequence"></textarea>   
		        </div>
		        <div class="col-md-2 col-sm-2 padd-0">
		        	<button class="default-btn btn-bg-two border-radius-50" type="submit" style="width:100%;padding: 5px;" onclick="return protein_click()">BLAST</button><hr>
					<p style="font-size:9px">An example protein sequence of  "<a style="text-decoration: none; cursor: pointer;margin-left:0px" onclick="document.getElementById('searchName1').value='EAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSRPRQVDHLRSGVQDSLANIAKSHLY'"><font color="#dcb200">TP53</font>".</a></p>
		        </div>
		        </div> 
		        </nav>
				</form>
		
		<hr>
							<iframe frameborder=0 width="100%" height="650px" 
							name="exmdb_blast" id="exmdb_blast"
							src="http://www.bio-server.cn/Server_blast/blast_res_new.jsp"
							scrolling="no" onload="load()"></iframe></div>
</div>



</div>
<footer class="footer-area footer-bg">
  <div class="container">
    <div class="copy-right-area">
      <div class="copy-right-text">
        <p>2022,CopyRight © HMU. College of Bioinformatics Science and Technology, Harbin, China.</p>
      </div>
    </div>
  </div>
</footer>
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

<!-- select -->
<script src="bootstrap_select/bootstrap-select.js"></script>


  <script>
    wow = new WOW(
      {
        animateClass: 'animated',
        offset:100
      }
    );
    wow.init();
  </script>

<script>
function load(){
	$("#graph_qi").fadeOut(1);
}

function gene_click(){
	if(document.getElementById("searchName").value==''){
	alert("Please input a sequence at least!");
	return false;
	}
	$("#graph_qi").fadeIn(1);
}

function protein_click(){
	if(document.getElementById("searchName1").value==''){
	alert("Please input a sequence at least!");
	return false;
	}
	$("#graph_qi").fadeIn(1);
}
function show_qi(){
	
}
</script>



</body>
</html>