<%@ page language="java" import="java.sql.*,wp.base.*" contentType="text/html; charset=UTF-8"
    pageEncoding="UTF-8"%>
<!DOCTYPE>
<head>
<title>ExMdb : browse</title>
<meta charset="utf-8">
<link rel="stylesheet" href="assets/css/bootstrap.min.css">
<link rel="stylesheet" type="text/css" href="select/demo.css"/>
<link rel="stylesheet" type="text/css" href="select/style-adsila.css" />
<link rel="stylesheet" href="select/selectpage_browse.css" type="text/css">


<link rel="stylesheet" type="text/css" href="tagtree/css/tagTree.css" />
<link rel="stylesheet" type="text/css" href="tagtree/css/font-awesome.min.css">
<script type="text/javascript" src="tagtree/js/jquery.min.js"></script>
<script type="text/javascript" src="tagtree/js/tagTree.js"></script>

<link rel="stylesheet" href="assets/css/style.min.css">

<script type="text/javascript">
<%
String search_table = dbhello.exmdb_left_list_exp();
String search_exp_list = dbhello.exmdb_left_exp();
String search_biomarker_list = dbhello.exmdb_left_biomarker();
// System.out.println(search_table2);
%>

$(function(){

    var data =[
    	{
	    	name:"Cancer Biomarkers",
	    	value:"all",
	    	children:[
	    		{name:"	Drug resistance [items = 904]",value:"	drug	",id:"marker",count:"535",children:[]},
	    		{name:"	Circulating [items = 878]",value:"	circulating	",id:"marker",count:"3449",children:[]},
	    		{name:"	Survival [items = 4263]",value:"	survival	",id:"marker",count:"4274",children:[]},
	    		{name:"	Immune [items = 1305]",value:"	immune	",id:"marker",count:"120",children:[]},
	    		{name:"	Metastasis [items = 3887]",value:"	metastasis	",id:"marker",count:"3278",children:[]},
	    		{name:"	Recurrence [items = 516]",value:"	recurrence	",id:"marker",count:"292",children:[]},
	    		{name:"	Cell Growth [items = 3854]",value:"	cellgrowth	",id:"marker",count:"2711",children:[]},
	    		{name:"	Epithelial Mesenchymal Transition [items = 3272]",value:"	emt	",id:"marker",count:"3038",children:[]},
	    		{name:"	Apoptosis [items = 3532]",value:"	apoptosis	",id:"marker",count:"2829",children:[]},
	    		{name:"	Cell Autophagy [items = 297]",value:"	autophagy	",id:"marker",count:"163",children:[]}
	    	]
	    },	    
	    {
	    	name:"Experimental Diseases",
	    	value:"all",
	    	children:[
				<%=search_table%>
	    	]
	    },
	    {
	    	name:"High-Throughput Bulk Datasets",
	    	value:"all",
	    	children:[
				{name:"GSE106817 [BRCA,C=115,N=2759]",value:"GSE106817_BRCA",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [COAD,C=115,N=2759]",value:"GSE106817_CRC",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [ESCA,C=88,N=2759]",value:"GSE106817_ESCA",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [STAD,C=115,N=2759]",value:"GSE106817_GC",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [HCC,C=115,N=2759]",value:"GSE106817_HCC",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [LUAD,C=115,N=2759]",value:"GSE106817_Lung",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [OV,C=320,N=2759]",value:"GSE106817_Ov",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [PAAD,C=115,N=2759]",value:"GSE106817_PAAD",id:"GSE_iframe",children:[]},
				{name:"GSE106817 [SARC,C=115,N=2759]",value:"GSE106817_SARC",id:"GSE_iframe",children:[]},
				{name:"GSE110271 [Multiple Myeloma,C=9,N=2]",value:"GSE110271",id:"GSE_iframe",children:[]},
				{name:"GSE112496 [ESCA,C=5,N=5]",value:"GSE112496",id:"GSE_iframe",children:[]},
				{name:"GSE112840 [ESCA,C=52,N=52]",value:"GSE112840",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [BRCA,C=40,N=100]",value:"GSE113486_BC",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [CHOL,C=40,N=100]",value:"GSE113486_Biliary",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [BLCA,C=392,N=100]",value:"GSE113486_Bladder",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [COAD,C=40,N=100]",value:"GSE113486_CRC",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [ESCA,C=40,N=100]",value:"GSE113486_Esophageal",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [STAD,C=40,N=100]",value:"GSE113486_GC",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [GBM,C=40,N=100]",value:"GSE113486_Glioma",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [HCC,C=40,N=100]",value:"GSE113486_Hepatocellular",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [LUAD,C=40,N=100]",value:"GSE113486_Lung",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [OV,C=40,N=100]",value:"GSE113486_Ov",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [PAAD,C=40,N=100]",value:"GSE113486_Pancreatic",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [PRAD,C=40,N=100]",value:"GSE113486_Prostate",id:"GSE_iframe",children:[]},
				{name:"GSE113486 [SARC,C=40,N=100]",value:"GSE113486_Sarcoma",id:"GSE_iframe",children:[]},
				{name:"GSE118782 [BRCA,C=30,N=10]",value:"GSE118782",id:"GSE_iframe",children:[]},
				
				{name:"GSE104926 [Esophagitis,C=6,N=6]",value:"GSE104926_esophagitis",id:"GSE_iframe",children:[]},
				{name:"GSE104926 [ESCA,C=6,N=6]",value:"GSE104926_ESCA",id:"GSE_iframe",children:[]},
				{name:"GSE106804 [GBM,C=8,N=6]",value:"GSE106804",id:"GSE_iframe",children:[]},
				
				{name:"GSE122488 [GBM,C=22,N=16]",value:"GSE122488",id:"GSE_iframe",children:[]},
				{name:"GSE122497 [ESCA,C=566,N=4965]",value:"GSE122497",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [BSTS,C=892,N=275]",value:"GSE124158_Bone",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [COAD,C=30,N=275]",value:"GSE124158_CRC",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [ESCA,C=30,N=275]",value:"GSE124158_ESCA",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [STAD,C=30,N=275]",value:"GSE124158_GC",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [GBM,C=30,N=275]",value:"GSE124158_Glioma",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [HCC,C=30,N=275]",value:"GSE124158_HCC",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [PAAD,C=30,N=275]",value:"GSE124158_Pancreatic",id:"GSE_iframe",children:[]},
				{name:"GSE124158 [LUAD,C=30,N=275]",value:"GSE124158_Lung",id:"GSE_iframe",children:[]},
				{name:"GSE125442 [KIRC,C=10,N=10]",value:"GSE125442",id:"GSE_iframe",children:[]},
				{name:"GSE126399 [STAD,C=12,N=10]",value:"GSE126399",id:"GSE_iframe",children:[]},
				{name:"GSE130512 [Tyriod Cancer,C=16,N=8]",value:"GSE130512",id:"GSE_iframe",children:[]},
				{name:"GSE130654 [STAD,C=36,N=12]",value:"GSE130654",id:"GSE_iframe",children:[]},
				{name:"GSE133684 [PAAD,C=284,N=117]",value:"GSE133684",id:"GSE_iframe",children:[]},
				{name:"GSE136321 [PRAD,C=24,N=24]",value:"GSE136321",id:"GSE_iframe",children:[]},
				{name:"GSE137140 [LUAD,C=2178,N=1746]",value:"GSE137140",id:"GSE_iframe",children:[]},
				{name:"GSE138740 [PRAD,C=146,N=89]",value:"GSE138740",id:"GSE_iframe",children:[]},
				{name:"GSE139031 [GBM,C=423,N=157]",value:"GSE139031",id:"GSE_iframe",children:[]},
				{name:"GSE144521 [CHOL,C=35,N=61]",value:"GSE144521",id:"GSE_iframe",children:[]},
				{name:"GSE149200 [COAD,C=5,N=5]",value:"GSE149200",id:"GSE_iframe",children:[]},
				{name:"GSE159177 [PRAD,C=278,N=187]",value:"GSE159177",id:"GSE_iframe",children:[]},
				{name:"GSE164174 [STAD,C=1423,N=1417]",value:"GSE164174_GC",id:"GSE_iframe",children:[]},
				{name:"GSE164174 [ESCA,C=1423,N=1417]",value:"GSE164174_ESCA",id:"GSE_iframe",children:[]},
				{name:"GSE164174 [COAD,C=1423,N=1417]",value:"GSE164174_CRC",id:"GSE_iframe",children:[]},
				{name:"GSE37472 [Oral Cancer,C=50,N=26]",value:"GSE37472",id:"GSE_iframe",children:[]},
				{name:"GSE41922 [BRCA,C=32,N=22]",value:"GSE41922",id:"GSE_iframe",children:[]},
				{name:"GSE63108 [ESCA,C=28,N=19]",value:"GSE63108",id:"GSE_iframe",children:[]},
				{name:"GSE67075 [COAD,C=18,N=30]",value:"GSE67075",id:"GSE_iframe",children:[]},
				{name:"GSE73002 [BRCA,C=1280,N=2686]",value:"GSE73002",id:"GSE_iframe",children:[]},
				{name:"GSE76449 [OV,C=24,N=4]",value:"GSE76449",id:"GSE_iframe",children:[]},
				{name:"GSE85589 [COAD,C=5,N=19]",value:"GSE85589_CC",id:"GSE_iframe",children:[]},
				{name:"GSE85589 [STAD,C=7,N=19]",value:"GSE85589_SC",id:"GSE_iframe",children:[]},
				{name:"GSE85589 [PAAD,C=88,N=19]",value:"GSE85589_PC",id:"GSE_iframe",children:[]},
				{name:"GSE85589 [CHOL,C=101,N=19]",value:"GSE85589_ICC",id:"GSE_iframe",children:[]},
	    	]
	    },

	    
	    {
	    	name:"Single-Cell Datasets",
	    	value:"all",
	    	children:[
	    		{name:"BLCA_GSE145281_aPDL1 [ N = 10 ] ",value:"BLCA_GSE145281_aPDL1",id:"sc_data",children:[]},
	    		{name:"CLL_GSE111014 [ N = 4 ] ",value:"CLL_GSE111014",id:"sc_data",children:[]},
	    		{name:"KIRC_GSE145281_aPDL1 [ N = 4 ] ",value:"KIRC_GSE145281_aPDL1",id:"sc_data",children:[]},
	    		{name:"MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] ",value:"MCC_GSE117988_aPD1aCTLA4",id:"sc_data",children:[]},
	    		{name:"MCC_GSE118056_aPDL1 [ N = 1 ] ",value:"MCC_GSE118056_aPDL1",id:"sc_data",children:[]},
	    		{name:"NSCLC_GSE127471 [ N = 1 ] ",value:"NSCLC_GSE127471",id:"sc_data",children:[]},
	    		{name:"PBMC_30K_10X [ N = 1 ] ",value:"PBMC_30K_10X",id:"sc_data",children:[]},
	    		{name:"PBMC_8K_10X [ N = 1 ] ",value:"PBMC_8K_10X",id:"sc_data",children:[]},
	    		{name:"PC_GSE67980 [ N = 13 ] ",value:"PC_GSE67980",id:"sc_data",children:[]},
	    		{name:"MEL_GSE157743 [ N = 22 ] ",value:"MEL_GSE157743",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE168212 [ N = 1 ] ",value:"COVID19_GSE168212",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE166992 [ N = 9 ] ",value:"COVID19_GSE166992",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE154567_severe [ N = 3 ] ",value:"COVID19_GSE154567_severe",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE154567_moderate [ N = 3 ] ",value:"COVID19_GSE154567_moderate",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE154567_recovery [ N = 3 ] ",value:"COVID19_GSE154567_recovery",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE167118_BALF_moderate [ N = 2 ] ",value:"COVID19_GSE167118_BALF_moderate",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE167118_BALF_severe [ N = 6 ] ",value:"COVID19_GSE167118_BALF_severe",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE167118_blood_moderate [ N = 1 ] ",value:"COVID19_GSE167118_blood_moderate",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE167118_blood_severe [ N = 6 ] ",value:"COVID19_GSE167118_blood_severe",id:"sc_data",children:[]},
	    		{name:"BP_GSE167118_BALF_moderate [ N = 3 ] ",value:"BP_GSE167118_BALF_moderate",id:"sc_data",children:[]},
	    		{name:"BP_GSE167118_BALF_no [ N = 1 ] ",value:"BP_GSE167118_BALF_no",id:"sc_data",children:[]},
	    		{name:"BP_GSE167118_blood_moderate [ N = 3 ] ",value:"BP_GSE167118_blood_moderate",id:"sc_data",children:[]},
	    		{name:"BP_GSE167118_blood_no [ N = 1 ] ",value:"BP_GSE167118_blood_no",id:"sc_data",children:[]},
	    		{name:"COVID19_GSE169503 [ N = 4 ] ",value:"COVID19_GSE169503",id:"sc_data",children:[]},
	    		{name:"SLE_GSE142016 [ N = 3 ] ",value:"SLE_GSE142016",id:"sc_data",children:[]},
	    		{name:"LIHC_GSE107747 [ N = 2 ] ",value:"LIHC_GSE107747",id:"sc_data",children:[]},
	    		{name:"Lymphoma_GSE124899 [ N = 2 ] ",value:"Lymphoma_GSE124899",id:"sc_data",children:[]},
	    		{name:"IPEX_GSE167976 [ N = 4 ] ",value:"IPEX_GSE167976",id:"sc_data",children:[]},
	    	]
	    },
    ];

    $("#test").tagTree({
    	id:"",
    	data:data,
    	fold:true,
    	multiple:false,
    	check:function(val){
    		console.log('chekc:'+val);
    		console.log($(this).tagTreeValues());
    	},
    	done:function(){
    		console.log('tagTree is ok!');
    	}
    });
});
</script>

<script type="text/javascript" src="select/demo2.js" ></script>
<script type="text/javascript" src="select/data/cell_out.js" ></script>
<script type="text/javascript" src="select/selectpage.min.js" ></script>
<script type="text/javascript">
	$(function(){
		var tag_data = [
			<%=search_exp_list%>
		];

		var tag_data2 = [
			<%=search_biomarker_list%>
		];


		$('#searchname').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data,
			lang: 'en',
			noResultClean : true
		});
		
		$('#Cell_location_searchname').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data_cell,
			lang: 'en',
			noResultClean : true
		});
		
		SyntaxHighlighter.all();
	});
</script>

<style>
::-webkit-scrollbar-track
{
	-webkit-box-shadow: inset 0 0 6px rgba(0,0,0,0.1);
	background-color: #F5F5F5;
	border-radius: 5px;
}

::-webkit-scrollbar
{
	width: 15px;
	background-color: #F5F5F5;
}

::-webkit-scrollbar-thumb
{
	border-radius: 10px;
	width: 15px;
	height: 5px;
	background-color: #E4665C;
}

div.sp_clear_btn{
	padding: 10px 0 0 0
}

div.sp_clear_btn i{
	font-size:14px
}

.sp_container{
	width:80% !important;
}


.tag_input{
	margin-left: 7px;
	font-size: 20px;
	float:left;
	margin-right:10px
}

.tag_input:before{
	content:" ";
    border-left:1px dashed #795548;
    bottom:50px;
    height:100%;
    top:-4;
    width:1px
}

#search_con
{
	margin-top: 2%;
	margin-bottom: 5px;
}

.label_tag{
    position: relative;
    z-index: 1;
    background-color: white;
    margin:5px 5px 5px 20px;
    /* height: 28px; */
    display: inline-block;
    padding: 0 8 0 5;
    font-weight: 700;
    font-size: large;
    line-height: 28px;
    color: black;
    transition: all 0.3s;
	-moz-transition: all 0.3s; /* Firefox 4 */
	-webkit-transition: all 0.3s; /* Safari 和 Chrome */
	-o-transition: all 0.3s; /* Opera */
}

.label_tag:before{
	content:" ";
	position: absolute;
    z-index: -1;
    left: -15px;
    border-top: 14px solid transparent;
    border-right: 14px solid #80b3ff;
    border-bottom: 14px solid transparent;
}


.label_tag:hover::before{
	border-right: 14px solid #80b3ff;
}

.label_tag:hover{
	background: #80b3ff;
    color: white;
    display: inline-block;
    cursor: pointer;
}


</style>


</head>

<body>
	
	<div id="test">
	<p class="label_tag" onclick="dieorlive()">Quick Search <i class="fa fa-check" id="tag_icon" style="display: none;"></i></p> <div class="node-count">117740</div>
	<div class="row" id="search_con" style="display:none">
	<div class="col-lg-10 col-md-10"><p class="tag_input">-- </p><input onchange="changeurl_exp(this.value)" type="text" name="searchname" id="searchname" class="form-control" value="MALAT1" placeholder=" XIST / MALAT1 / ..." style="font-size:14px;font-family: 'Montserrat';">
	</div>
	</div>
	
	<p> </p>
	<p class="label_tag" onclick="dieorlive_cell_location()" style="margin-bottom: 12px;">Cellular Localization <i class="fa fa-check" id="tag_icon_cell_location" style="display: none;"></i></p> <div class="node-count">56098</div>
	<div class="row" id="search_con_cell_location" style="display:none">
	<div class="col-lg-10 col-md-10"><p class="tag_input">-- </p><input onchange="changeurl_cell_location(this.value)" type="text" name="Cell_location_searchname" id="Cell_location_searchname" class="form-control" value="MALAT1" placeholder=" XIST / MALAT1 / ..." style="font-size:14px;font-family: 'Montserrat';">
	</div>
	</div>
	
	
	
<!-- 	<h3>Biomarker Genes</h3> -->
<!-- 	<input onchange="changeurl_bio(this.value)" type="text" name=searchname id="searchname2" class="form-control" placeholder=" XIST / MALAT1 / ..." style="font-size:30px;padding:30px 15px;font-family: 'Montserrat';"> -->
	
	</div>
	

	<div class="browse_right" id="browse_right">
		<iframe frameborder=0 name="right_iframe" id="right_iframe" style="width: 100%; height: 900px;margin: 0 auto; position: relative;" src="exmdb_browse_panel.jsp?searchname=MALAT1";> </iframe>
	</div>
</body>


<script>
function changeurl_exp(val)
{
	document.getElementById("right_iframe").src="exmdb_browse_panel.jsp?searchname="+val;
}

function changeurl_cell_location(val)
{
	document.getElementById("right_iframe").src="cellmap.jsp?searchname="+val;
}

function dieorlive()
{
	
	$("#search_con").slideToggle();
	$("#tag_icon").slideToggle();
	
}


function dieorlive_cell_location()
{
	$("#search_con_cell_location").slideToggle();
	$("#tag_icon_cell_location").slideToggle();
	
}
</script>

</html>