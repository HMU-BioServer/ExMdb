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
<title>Techex - Technology Services HTML Template</title>

<%
String searchname = "XIST";  

HashMap hm = new HashMap();
hm.put("predict",0); 
hm.put("biomarker",0); 
hm.put("celllocation",0); 
hm.put("exp",0); 

//该条目涉及了几个疾病
HashMap hm_dis = new HashMap();
hm_dis.put("predict",0); 
hm_dis.put("biomarker",0); 
hm_dis.put("celllocation",0); 
hm_dis.put("exp",0); 

//该条目涉及了几个样本
HashMap hm_sample = new HashMap();
hm_sample.put("predict",0); 
hm_sample.put("biomarker",0); 
hm_sample.put("celllocation",0); 
hm_sample.put("exp",0); 

//该条目涉及了几个样本
HashMap hm_data = new HashMap();
hm_data.put("predict",0); 
hm_data.put("biomarker",0); 
hm_data.put("celllocation",0); 
hm_data.put("exp",0); 

if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}


String Zero_stop = "";

String qi="";
qi = dbhello.quickString(searchname);
// System.out.println(qi);

if(!qi.equals("")){
	for (String retval: qi.split(",")){
		
		String[] qiyue = retval.split("_");
		hm.put(qiyue[0],qiyue[1]); 
		hm_dis.put(qiyue[0],qiyue[2]); 
		hm_sample.put(qiyue[0],qiyue[3]); 
		hm_data.put(qiyue[0],qiyue[4]); 
	}	
}

String trans_searchname = dbhello.exmdb_quick_trans(searchname);
String trans_count = "0";

if (!trans_searchname.equals("")&&trans_searchname!=null&&!trans_searchname.equals("NA")){
	trans_count = dbhello.exmdb_quick_trans_count(trans_searchname);
}else{}


String sc_res = dbhello.exmdb_get_quick_count(searchname);

HashMap sc_data = new HashMap();
sc_data.put("diseasename1",0); 
sc_data.put("datasetname",0); 
sc_data.put("allcount",0); 


if(!sc_res.equals("")){
	String[] retval = sc_res.split("_");
	sc_data.put("diseasename1",retval[0]); 
	sc_data.put("datasetname",retval[1]); 
	sc_data.put("allcount",retval[2]); 
}

%>

<style>
body{
	overflow: hidden;
}
}
.search_item{
	font-size:20px;
}

.table_contact{
	padding: 1rem;
	font-size: 18px;
}

.sub_num{
	font-weight: 900;
	color: #265cd4;
}

.sub_num:hover{
	color: #ffc221;
	-webkit-transition:all .3s;
	transition:all .3s;
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








<div class="terms-conditions-area" style="padding-top:10px">
  <div class="container" style="max-width: 1440px;">
  <h4>Records in ExMdb : <span style="color: #1296db;"><%=searchname%></span> [ <span style="color: #a4e0ff;"><%=trans_searchname%></span> ]</h4>
  <hr>
	<div class="row">
	      <div class="col-lg-3 col-md-6 wow fadeInUp animated" data-wow-delay="0.3s">
	        <div class="blog-card">
	          <div class="blog-img"> <a href="exmdb_quick_high_data.jsp?searchname=<%=searchname%>" target="_blank"> <img src="assets/images/quick_pic/Highput_4.png" alt="Blog Images"> </a>
	          </div>
	          <div class="content">
	            <h3 class="center_text"> <a target="_blank" href="exmdb_quick_high_data.jsp?searchname=<%=searchname%>">High-throughput entries</a> </h3>
<!-- 	            <img style="width:86.5%" src="assets/images/icon/renti.png" alt="Blog Images"> -->

	            <p class="search_item">
	            
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#30457b" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#f5b90f" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class=""></path></svg>
	             Entries Counts : <a target="_blank" href="exmdb_quick_high_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=Integer.valueOf(trans_count)+Integer.parseInt(hm.get("predict").toString())%></span></a></p>
	            
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#f5b90f" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#30457b" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class="selected"></path></svg>	             
	              No. of Diseases : <a target="_blank" href="exmdb_quick_high_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=hm_dis.get("predict") %></span></a></p>

	             <p class="search_item">
	             <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#30457b" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#f5b90f" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class=""></path></svg>
	             No. of DataSet : <a target="_blank" href="exmdb_quick_high_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=Integer.valueOf(trans_count)+Integer.parseInt(hm.get("predict").toString())%></span></a></p>
	           </div>
	        </div>
	      </div>
	      
	      <div class="col-lg-3 col-md-6 wow fadeInUp animated" data-wow-delay="0.6s">
	        <div class="blog-card">
	          <div class="blog-img"> <a href="exmdb_quick_exp_data.jsp?searchname=<%=searchname%>" target="_blank"> <img src="assets/images/quick_pic/exp_2.png" alt="Blog Images"> </a>

	          </div>
	          <div class="content">
	            <h3 class="center_text"> <a target="_blank" href="exmdb_quick_exp_data.jsp?searchname=<%=searchname%>">Experimental validation</a> </h3>
	            
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#f5b90f" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#30457b" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class="selected"></path></svg>	             
	             All Counts : <a target="_blank" href="exmdb_quick_exp_data.jsp?searchname=<%=searchname%>"> <span class="sub_num"><%=hm.get("exp") %></span></a></p>
	            
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#30457b" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#f5b90f" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class=""></path></svg>
	             No. of Disease : <a target="_blank" href="exmdb_quick_exp_data.jsp?searchname=<%=searchname%>"> <span class="sub_num"><%=hm_dis.get("exp") %></span></a></p>
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#f5b90f" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#30457b" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class="selected"></path></svg>	             
	             No. of Species : <a target="_blank" href="exmdb_quick_exp_data.jsp?searchname=<%=searchname%>"> <span class="sub_num"><%=hm_data.get("exp") %></span></a></p>
	          </div>
	        </div>
	      </div>
	      
	      <div class="col-lg-3 col-md-6 wow fadeInUp animated" data-wow-delay="0.9s">
	        <div class="blog-card">
	          <div class="blog-img"> <a href="exmdb_quick_marker_data.jsp?searchname=<%=searchname%>" target="_blank"> <img src="assets/images/quick_pic/biomarker_2.png" alt="Blog Images"> </a>

	          </div>
	          <div class="content">

	            <h3 class="center_text"> <a target="_blank" href="exmdb_quick_marker_data.jsp?searchname=<%=searchname%>">Cancer biomarkers</a> </h3>
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#30457b" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#f5b90f" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class=""></path></svg>
	            Biomarker Counts: <a target="_blank" href="exmdb_quick_marker_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=hm.get("biomarker") %></span></a></p>
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#f5b90f" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#30457b" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class="selected"></path></svg>	            
	             No. of Disease : <a target="_blank" href="exmdb_quick_marker_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=hm_dis.get("biomarker") %></span></a></p>
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#30457b" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#f5b90f" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class=""></path></svg>
	             No. of Studies : <a target="_blank" href="exmdb_quick_marker_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=hm_data.get("biomarker") %></span></a></p>
	          </div>
	        </div>
	      </div>
	      
	      <div class="col-lg-3 col-md-6 wow fadeInUp animated" data-wow-delay="1.2s">
	        <div class="blog-card">
	          <div class="blog-img"> <a href="exmdb_quick_sc_data.jsp?searchname=<%=searchname%>" target="_blank"> <img src="assets/images/quick_pic/sc_res_4.png" alt="Blog Images"> </a>
	          </div>
	          <div class="content">

	            <h3 class="center_text"> <a target="_blank" href="exmdb_quick_sc_data.jsp?searchname=<%=searchname%>">Single Cell Analysis</a> </h3>
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#f5b90f" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#30457b" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class="selected"></path></svg>	             
	            Disease Counts: <a target="_blank" href="exmdb_quick_sc_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=sc_data.get("diseasename1") %></span></a></p>
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#30457b" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#f5b90f" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class=""></path></svg>
	            No. of Data Set : <a target="_blank" href="exmdb_quick_sc_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=sc_data.get("datasetname") %></span></a></p>
	            <p class="search_item">
	            <svg style="margin-top: -2%;" t="1652689225787" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2644" data-spm-anchor-id="a313x.7781069.0.i11" width="30" height="30"><path d="M259.981 512.012l252.026-252.021 251.997 251.997-252.021 252.021-252.002-251.997z m0 0" fill="#f5b90f" p-id="2645" data-spm-anchor-id="a313x.7781069.0.i0" class=""></path><path d="M511.997 64L64 511.997 511.997 960 960 511.997 511.997 64z m0 806L154 511.997 511.997 154 870 511.997 511.997 870z m0 0" fill="#30457b" p-id="2646" data-spm-anchor-id="a313x.7781069.0.i6" class="selected"></path></svg>	             
	            All Count : <a target="_blank" href="exmdb_quick_sc_data.jsp?searchname=<%=searchname%>"><span class="sub_num"><%=sc_data.get("allcount") %></span></a></p>
	           </div>
	        </div>
	      </div>
	</div>

	<div class="contact-wrapper wow fadeInUp animated" data-wow-delay="1.6s" style="padding:20px;">
    	<h4 class="home_title" style="color:#424F60">External Annotation for <span style="color:#F9B189;text-transform:uppercase"><%=searchname%></span></h4><hr>
    	
    	<div class="row">
    		<div class="col-xl-6 col-lg-6 col-md-6 col-sm-6" style="border-right: 2px solid #cacaca;">
    		
		    	<div class="table_contact">
				<div class="row" style="height:100%">
					<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="http://www.noncode.org/gene_trans_search.php?keyword=<%=searchname%>"><img class='table_img' style='width:200px;' src="img/links/noncode.png"></a></div>
					<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">An integrated knowledge database dedicated to ncRNAs, especially lncRNAs.</div>
				</div>
				</div>
				
				<div class="table_contact">
				<div class="row" style="height:100%">
					<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box" ><a href="https://pubmed.ncbi.nlm.nih.gov/?term=<%=searchname%>"><img class='table_img' style='width:190px;' src="img/links/pubmed.png"></a></div>
					<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">PubMed® comprises more than 32 million citations for biomedical literature from MEDLINE, life science journals, and online books.</div>
				</div>
				</div>
				
				
				<div class="table_contact">
				<div class="row" style="height:100%">
					<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="https://www.ncbi.nlm.nih.gov/genome/?term=<%=searchname%>"><img class='table_img' style='width:190px;' src="img/links/NCBI.png"></a></div>
					<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">National Center for Biotechnology Information.</div>
				</div>
				</div>
    		
    		</div>
    		<div class="col-xl-6 col-lg-6 col-md-6 col-sm-6">
    				<div class="table_contact">
					<div class="row" style="height:100%">
						<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="https://www.gencodegenes.org/"><img class='table_img' style='width:190px;' src="img/links/Gencode.png"></a> 	</div>
						<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">GENCODE are supporting the annotation of non-canonical human ORFs predicted by Ribo-seq data.</div>
					</div>
					</div>
			
					<div class="table_contact">
					<div class="row" style="height:100%">
						<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box" ><a href="https://www.genecards.org/lookup/<%=searchname%>"><img class='table_img' style='width:190px;' src="img/links/genecard.png"></a> 	</div>
						<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box" >GeneCards enables researchers to effectively navigate and inter-relate the wide universe of human genes, diseases etc.</div>
					</div>
					</div>
			
					<div class="table_contact">
					<div class="row" style="height:100%">
						<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="https://omim.org/search?search=<%=searchname%>"><img class='table_img' style='width:190px;' src="img/links/OMIM.png"></a></div>
						<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">A Knowledgebase of Human Genes and Genetic Phenotypes</div>
					</div>
					</div>
    		</div>
    	</div>
    </div>

  </div>
</div>


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