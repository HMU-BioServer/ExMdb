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


<title>ExMdb : Function enrichment tools</title>



<style>
.search_item{
	font-size:20px;
}

.imm_select{
	width:100%;
}


.theme-btn {
	margin-top:30%;
	width:50%;
	float:left;
	padding:4%;
	font-size: 15px;
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

div.sp_clear_btn {
    font-size: 25px;
}

ul.sp_element_box{
	padding: 5px 0 0 3px;
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
  
</div>

<div class="inner-banner">
  <div class="container">
    <div class="inner-title text-center">
      <h3>Function enrichment tools</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Function enrichment tools</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
<div class="terms-conditions-area" style="padding:40px 0px">
  <div class="container wow fadeInUp animated" style="max-width: 1440px;">

		<div class="contact-header">
			<h1>Function Enrichment Panel</h1>
		</div>	
		<div class="contact-wrapper">
	    <form action="exmdb_function_mRNA_in.jsp" target="exmdb_function"  method="post" role="form" id="imm_form" name="imm_form">
	    <p class="sub_title_function"><img src="img/iconfont/function.png" style="margin-top: -7px;margin-right: 7px;">ExMdb-Function</p>
	    <p class="sub_text_function">Identify dysregulated functions of mRNA based on Gene Ontology and biological pathways.</p>
            <div class="row">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-10 col-sm-10 col-xs-10">
							<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Genes</span>
									<span class="menu__item-label">(Enter keywords)</span>
								</a>
								<input type="text" name="searchname" id="searchname" class="form-control" value="TP53,CDK6,EGFR,VEGFA" placeholder=" TP53 / BRCA1 / ..."  style="font-size:15px;padding:20px 15px">
							</nav>
            </div></div>
            <div class="row">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-6 col-sm-6 col-xs-6">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">No. of top enrichmented functions</span>
									<span class="menu__item-label">(Display Counts)</span>
								</a>
								<select class="language-list-item imm_select" name="top" id="top">
					              <option value="5">5</option>
					              <option value="10" selected="selected">10</option>
					              <option value="20">20</option>
					            </select>
							</nav>

            </div>
            <div class="col-md-2 col-sm-2 col-xs-2">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Color</span>
									<span class="menu__item-label">Upper / Lower</span>
								</a>
								<input id="top_color" name="top_color" type="color" value="#cca300" style="height:40px;width:80px">  <input name="bottom_color" id="bottom_color" type="color" value="#002db3" style="height:40px;width:80px">
							</nav>

            </div>
            
            <div class="col-md-2 col-sm-2 col-xs-2">
            	<button class="default-btn btn-bg-two border-radius-50 theme-btn" type="button" onclick="return gene_random()">EXAMPLE</button>
            	<button class="default-btn btn-bg-two border-radius-50 theme-btn" id="submit_botton" onclick="show_qi()" type="submit">SUBMIT</button>
            </div>
            </div>
        </form>
        <br>
        <hr>
						<iframe frameborder=0 width="100%" height="1300px"
							name="exmdb_function" id="exmdb_function"
							src="exmdb_function_mRNA_in.jsp"
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


<script type="text/javascript" src="select/demo2.js" ></script>
<script type="text/javascript" src="select/jquery.min.js" ></script>
<script type="text/javascript" src="select/selectpage.min.js" ></script>
<script type="text/javascript" src="select/data/mRNA_list.js" ></script>
<script type="text/javascript">
	$(function(){
		$('#searchname').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data,
			lang: 'en',
			multiple : true,
			maxSelectLimit : 40,
			noResultClean : true
		});
		
	
		SyntaxHighlighter.all();
	});
</script>

<script language="javascript">
	function gene_random() {
//		var array = new Array('VEGFA', 'PTEN','IGF1R', 'CDK6', 'SMAD4','TP53','EGFR','NFKB1','ENSG00000112715', 'ENSG00000171862','ENSG00000140443','ENSG00000105810','ENSG00000141646');
		document.getElementById('searchname').value = "TP53,CDK6,EGFR,VEGFA";
		document.getElementById("submit_botton").click();
	}
	function reset_quick() {
		document.getElementById('searchname').value = "";
	}	
	
	
	function load(){
		$("#graph_qi").fadeOut(1);
	}

	function show_qi(){
		$("#graph_qi").fadeIn(1);
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