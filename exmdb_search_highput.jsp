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
<link rel="stylesheet" href="select/selectpage_cell.css" type="text/css">


<link rel="stylesheet" type="text/css" href="scroll_bar/css/styles.css">
<link rel="stylesheet" type="text/css" href="scroll_bar/css/jquery-ui.min.css">

<title>ExMdb : Search [High throughput]</title>

<%
String searchname = "XIST";  



if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}



%>

<style>
.search_item{
	font-size:20px;
}

div.sp_clear_btn {
    font-size: 5px;
}

div.sp_clear_btn i {
    font-size: 22px;
}

.imm_select{
	width:100%;
}
.list{
	height: 120px;
    overflow-y: auto !important;
}

.nice-select{
	font-size:17px
}
}
</style>

</head>
<body>
<div class="preloader">
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
            <li class="nav-item"><a href="#" class="nav-link active">Search <i class='bx bx-caret-down'></i></a>
              <ul class="dropdown-menu">
                <li class="nav-item"><a href="exmdb_search_highput.jsp" class="nav-link" target="_blank">Search [ High throughput ]</a></li>
                <li class="nav-item"><a href="exmdb_search_exp.jsp" class="nav-link" target="_blank">Search [ Experimental validation ]</a></li>
                <li class="nav-item"><a href="exmdb_search_biomarker.jsp" class="nav-link" target="_blank">Search [ Cancer Biomarkers ]</a></li>
                <li class="nav-item"><a href="exmdb_search_scdataset.jsp" class="nav-link" target="_blank">Search [ Single-cell ]</a></li>
              </ul>
            </li>
            <li class="nav-item"><a href="#" class="nav-link">Tools <i class='bx bx-caret-down'></i></a>
              <ul class="dropdown-menu">
                <li class="nav-item"><a href="exmdb_cell_out.jsp" target="_blank" class="nav-link">ExMdb-Cellloaction </a></li>
                <li class="nav-item"><a href="exmdb_imm_out.jsp" target="_blank" class="nav-link">ExMdb-Immunity </a></li><li class="nav-item"><a href="#" class="nav-link">ExMdb-Survival Tools<i class='bx bx-caret-down'></i></a><ul class="dropdown-menu"><li class="nav-item"><a href="exmdb_survival_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">mRNA-Survival Analyze</a></li><li class="nav-item"><a href="exmdb_survival_miRNA_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">miRNA-Survival Analyze</a></li> <li class="nav-item"><a href="exmdb_survival_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">lncRNA-Survival Analyze</a></li></ul></li>
                <li class="nav-item"><a href="exmdb_net_out.jsp" target="_blank" class="nav-link">ExMdb-Network </a></li><li class="nav-item"><a href="exmdb_blast_out.jsp" target="_blank" class="nav-link">ExMdb-BLAST </a></li><li class="nav-item"><a href="exmdb_function_time.jsp" target="_blank" class="nav-link">ExMdb-Pseudotime pathway</a></li><li class="nav-item"><a href="exmdb_single_vis.jsp" target="_blank" class="nav-link">ExMdb-Single Cell Visualization</a></li>
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
  <div class="side-nav-responsive">
    <div class="container-max">
      <div class="dot-menu">
        <div class="circle-inner">
          <div class="in-circle circle-one"></div>
          <div class="in-circle circle-two"></div>
          <div class="in-circle circle-three"></div>
        </div>
      </div>
      <div class="container">
        <div class="side-nav-inner">
          <div class="side-nav justify-content-center align-items-center">
            <div class="side-nav-item nav-side">
              <div class="search-box"><i class='bx bx-search'></i></div>
              <div class="get-btn"><a href="contact.html" class="default-btn btn-bg-two border-radius-50">Get A Quote <i class='bx bx-chevron-right'></i></a></div>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</div>

<div class="inner-banner">
  <div class="container">
    <div class="inner-title text-center">
      <h3>Search [High throughput]</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Search [High throughput]</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
<div class="terms-conditions-area" style="padding:40px 0px">
<div class="container wow fadeInUp animated" style="max-width:1440px">

<div class="contact-header">
<h1>Search Panel</h1>
</div>
<div class="contact-wrapper">
<form id="gene" name="gene" method="post" target="_blank" action="exmdb_quick_high_data.jsp">
<div class="row">
<h4><i class="bx bx-star"></i> Quick Search</h4>
<hr>
	
	<div class="col-lg-1 col-md-1"></div>
	<div class="col-lg-6 col-md-6">
					<nav class="menu menu--adsila">
					<a class="menu__item" href="#"><span class="menu__item-name">Enter key words</span></a>	
					
						<input type="text" name="searchname" id="searchname" class="form-control" value="MALAT1" placeholder=" XIST / MALAT1 / ..." style="font-size:20px;">
					
					</nav>
	
	</div>
	<div class="col-lg-2 col-md-2">
					<button class="submit_my" type="button" onclick="return gene_random()"><span>Example</span></button>
	</div>
	<div class="col-lg-2 col-md-2">
					<button class="submit_my" type="submit" onclick="return gene_click()"><span>Submit</span></button>
	</div>
</div></form>

					
					
					<hr>
					<br>
					<form id="gene_ad" name="gene_ad" method="post" action="exmdb_search_table.jsp" target="_blank">
								<div class="row">
								<h4><i class="bx bx-search"></i> Advance Search</h4>
								<hr>
									<div class="col-lg-3 col-md-3" style="margin-left: 4%;">
									<nav class="menu menu--adsila">
										<a class="menu__item" href="#">
											<span class="menu__item-name">Enter Genes & Organ</span>
										</a>
										<input type="text" name="ADsearchname" id="ADsearchname" class="form-control" value="MIMAT0003236" placeholder=" XIST / MALAT1 / ..." style="font-size:20px;">
									<hr>
									<select class="language-list-item imm_select imm_select" name="organ" id="organ" >
										<option value="bone">Bone</option>
										<option value="bladder">Bladder</option>
										<option value="breast">Breast</option>
										<option value="cholangiocarcinoma ">Cholangiocarcinoma </option>
										<option value="colorectal">Colorectal</option>
										<option value="esophageal">Esophageal</option>
										<option value="gastric">Gastric</option>
										<option value="glioma">Glioma</option>
										<option value="liver">Liver</option>
										<option value="lung">Lung</option>
										<option value="oral ">Oral</option>
										<option value="ovarian">Ovarian</option>
										<option value="pancreatic">Pancreatic</option>
										<option value="prostate">Prostate</option>
										<option value="renal">Renal</option>
										<option value="tyriod">Tyriod</option>        
							         </select>	
							         </nav>
									</div>
									
									<div class="col-lg-4 col-md-4">	
									<nav class="menu menu--adsila">
										<a class="menu__item" href="#">
											<span class="menu__item-name">Drag the P-value bar</span>
										</a>					
									  <div class="cube">
									    <div class="a"></div>
									    <div class="b"></div>
									    <div class="c"></div>
									    <div class="d"></div>
									    <div id="slider-range-min"></div>
									  </div>
									 </nav>
									</div>
									<div class="col-lg-2 col-md-2">	
										<nav class="menu menu--adsila">
											<a class="menu__item" href="#">
												<span class="menu__item-name">P-value</span>
											</a>	
											<input type="text" id="amount" name="amount"/> 
										</nav>
									</div>
									
									<div class="col-lg-2 col-md-2">
										<button style="margin-top: 35%;font-size: 2rem;" class="submit_my" type="submit" onclick="return gene_click_AD()"><span>Submit</span></button>
									</div>
								</div>
					</form>
							
</div>
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

<script src="scroll_bar/js/jquery-2.1.1.min.js"></script>
<script src='scroll_bar/js/jquery-ui.min.js'></script>
	<script>
	$(function ($) {
	    $('#slider-range-min').slider({
	        range: 'min',
	        value: 0.5,
	        min: 0,
			step:0.01,
	        max: 1,
	        slide: function (event, ui) {
	            $('#amount').val(ui.value);
	            $('.a, .b, .c, .d').width(ui.value*100 + '%');
	        }
	    });
	    $('.ui-slider-handle').text(' ');
	    $('#amount').val($('#slider-range-min').slider('value'));
	});
</script>

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
<script type="text/javascript" src="select/data/highput.js" ></script>




<script type="text/javascript">
	$(function(){
		$('#searchname').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data,
			lang: 'en',
			noResultClean : true
		});
		
		$('#ADsearchname').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data,
			lang: 'en',
			noResultClean : true
		});
		
		SyntaxHighlighter.all();
	});
</script>


<script language="javascript">


	function gene_click(){
		if(document.getElementById("searchname").value==''){
		alert("please input a keyword at least!");
	return false;
	}}
	
	function gene_click_AD(){
		if(document.getElementById("ADsearchname").value==''){
		alert("please input a keyword at least!");
		return false;
	}}

	function gene_random() {
//		var array = new Array('VEGFA', 'PTEN','IGF1R', 'CDK6', 'SMAD4','TP53','EGFR','NFKB1','ENSG00000112715', 'ENSG00000171862','ENSG00000140443','ENSG00000105810','ENSG00000141646');
		var array = new Array('KRAS', 'ENSG00000133703','MALAT1','ENSG00000251562');

		var Ra_qi = Math.random();
		document.getElementById('searchname').value = array[Math.floor(Ra_qi* array.length)];
		document.getElementById('searchname_text').value = array[Math.floor(Ra_qi* array.length)];
	}
	function reset_quick() {
		document.getElementById('searchname').value = "";
	}	
</script>

<script>
    wow = new WOW(
      {
        animateClass: 'animated',
        offset:100
      }
    );
    wow.init();
</script>



</body>
</html>