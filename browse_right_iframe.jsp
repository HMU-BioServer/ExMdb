<%@ page language="java" contentType="text/html; charset=ISO-8859-1"
    pageEncoding="ISO-8859-1"%>
<!DOCTYPE html>
<html lang="zh">
<head>
	<meta charset="UTF-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1"> 
	<meta name="viewport" content="width=device-width, initial-scale=1.0">
	<title>Web</title>
	<link rel="stylesheet" type="text/css" href="right_file/css/normalize.css" /><!--CSS RESET-->
	<link rel="stylesheet" type="text/css" href="right_file/css/htmleaf-demo.css"><!--演示页面样式，使用时可以不引用-->
	<link rel="stylesheet" href="right_file/css/style.css">
<!-- 	<link href="http://cdn.bootcss.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet"> -->
	<link rel="stylesheet" type="text/css" href="assets/fonts/font-awesome.min.css" /><!--CSS RESET-->
	
	<style>
	
	.block_div{
	height:18%;
	width:80%;
	margin:auto;
	margin-bottom:10px;
	border-radius: 15px;
	padding:5px;
	}
	
	.block_left{
	width:30.5%;
	height:100%;
	float:left;
	border-radius:15px;
	padding-top: 8px;
	}
	
	.block_right{
	background:blue;
	width:69%;
	height:100%;
	float:right;
	border-radius:15px;
	}
	
	</style>
	
</head>
<body>
	<div class="wrapper">
           <div class="container">

               <div class="slideshow">

                    <div class="slideshow-left">

                        <div class="Lslide" style="background-image: url(assets/images/quick_pic/left_1.png);">
                            <div class="Lslide-content">
                                <h2>Data Statistics</h2>
                                <p>High throughput data / Experimental validation data</p>
								<img style="width:200px;background: #ffffff82;padding: 5%;border-radius:20px;border: 10px solid white;" src="img/browse_iframe/browse_1.png">

                            </div>
                        </div>
                        <div class="Lslide">
                            <div class="Lslide-content">
                                <h2>Disease entries</h2>
                                <p>Experimental confirmation entries</p>
								<img style="width:200px;background: #ffffff82;padding: 5%;border-radius:20px;border: 10px solid white;" src="img/browse_iframe/browse_2.png">
                            </div>
                        </div>
                        <div class="Lslide">
                            <div class="Lslide-content">
                                <h2>Disease Bodymap</h2>
                                <p>Diseases - Organs</p>
								<img style="width:200px;background: #ffffff82;padding: 5%;border-radius:20px;border: 10px solid white;" src="img/browse_iframe/browse_3.png">
                            </div>
                        </div>    
                        <div class="Lslide">
                            <div class="Lslide-content">
                                <h2>Analysis Tools</h2>
                                <p>Analysis and plotting</p>
								<img style="width:200px;background: #ffffff82;padding: 5%;border-radius:20px;border: 10px solid white;" src="img/browse_iframe/browse_4.png">
                            </div>
                        </div>   
                    </div>

                    <div class="slideshow-right">
                        <div class="Rslide">
                            <div id="container" style="height: 100%"></div>
                        </div>
                        <div class="Rslide">
                            <div id="container_bar" style="height: 100%"></div>
                        </div>     
                        <div class="Rslide">
                           <iframe frameborder=0  name="exmdb_imm" id="exmdb_imm" src="exmdb_body_table.html" style="width:1500px;height:100%" scrolling="no"></iframe>
                        </div>     
                        <div class="Rslide">
                           <div class="block_div" style="margin-top:10px">
                           	<div class="block_left"><img style="width:100%" src="assets/images/quick_pic/network.png"></div>
                           	<div class="block_right"></div>
                           </div>
                           <div class="block_div">
                            <div class="block_left"><img style="width:100%" src="assets/images/quick_pic/function.png"></div>
                           	<div class="block_right"></div>
                           </div>
                           <div class="block_div">
                           	<div class="block_left"><img style="width:100%" src="assets/images/quick_pic/celllocation.png"></div>
                           	<div class="block_right"></div>
                           </div>
                           <div class="block_div">
                           	<div class="block_left"><img style="width:100%" src="assets/images/quick_pic/imm_home.png"></div>
                           	<div class="block_right"></div>
                           </div>
                           <div class="block_div">
                           	<div class="block_left"><img style="width:100%" src="assets/images/quick_pic/sc_res_4.png"></div>
                           	<div class="block_right"></div>
                           </div>
                        </div>                                          
                    </div>    

                    
                    <div class="control">
                        <div class="oncontrol control-bottom">
                            <i class="fa fa-arrow-up" aria-hidden="true"></i>
                        </div>
                        <div class="oncontrol control-top">
                            <i class="fa fa-arrow-down" aria-hidden="true"></i>
                        </div>                          
                    </div>

               </div>

           </div>
       </div>
	<!-- <div class="htmleaf-container">
		<header class="htmleaf-header">
			<h1>js左右两侧分屏动画轮播图特效 <span>Slider boomerang effect</span></h1>
			<div class="htmleaf-links">
				<a class="htmleaf-icon icon-htmleaf-home-outline" href="http://www.htmleaf.com/" title="jQuery之家" target="_blank"><span> jQuery之家</span></a>
				<a class="htmleaf-icon icon-htmleaf-arrow-forward-outline" href="http://www.htmleaf.com/jQuery/Slideshow-Scroller/201805225134.html" title="返回下载页" target="_blank"><span> 返回下载页</span></a>
			</div>
		</header>
	</div> -->
	
	<script type="text/javascript">
		var Lslide      = document.querySelectorAll('.Lslide'),
		    Rslide      = document.querySelectorAll('.Rslide'),
		    control     = document.querySelectorAll('.oncontrol'),
		    slideHeight = document.querySelector('.wrapper').clientHeight,
		    color = ['#fdc97c', 'rgb(255 108 83)', 'rgb(81 195 255)','#263a9a'],
		    index = 0;


		function init() {
		    slideHeight = document.querySelector('.wrapper').clientHeight;
		    for (i = 0; i < Lslide.length; i++) {
		        Lslide[i].style.backgroundColor = color[i];
		        Lslide[i].style.top = -slideHeight * i + "px";
		        Rslide[i].style.top = slideHeight * i + "px";
		    }  
		}
		init()
		window.addEventListener('resize', function(){
		    init()
		});

		function moveToTop() {

		    index++;
		    for (el = 0; el < Lslide.length; el++) {
		        Lslide[el].style.top = parseInt(Lslide[el].style.top) + slideHeight + "px";
		        Rslide[el].style.top = parseInt(Rslide[el].style.top) - slideHeight + "px";
		    }

		    if (index > Lslide.length-1) {
		        index = 0;
		        for (el = 0; el < Lslide.length; el++) {
		            Lslide[el].style.top = -slideHeight * el + "px";
		            Rslide[el].style.top = slideHeight * el + "px";
		        }
		    }
		}

		function moveToBottom() {
		    index--;
		    for (el = 0; el < Lslide.length; el++) {
		        Lslide[el].style.top = parseInt(Lslide[el].style.top) - slideHeight + "px";
		        Rslide[el].style.top = parseInt(Rslide[el].style.top) + slideHeight + "px";
		        
		    }
		    if (index < 0) {
		        index = Rslide.length-1;
		        for (el = 0; el < Lslide.length; el++) {
		            Lslide[el].style.top = parseInt(Lslide[el].style.top) + slideHeight * Lslide.length + "px";
		            Rslide[el].style.top = parseInt(Rslide[el].style.top) - slideHeight * Rslide.length + "px";
		        }
		    }
		}

		function transition() {
		    for (t = 0; t < Lslide.length; t++) {
		        Lslide[t].style.transition = "all 0.8s";
		        Rslide[t].style.transition = "all 0.8s";
		    }
		}
		  

		for (t = 0; t < control.length; t++) {
		    control[t].addEventListener("click", function() {

		        if (this.classList.contains('control-top')) {
		            moveToTop()
		        } 
		        if (this.classList.contains('control-bottom')) {
		            moveToBottom()
		        }

		        transition()
		   
		    });
		}

		var wheeling;
		function mousemouve(e) {

		    clearTimeout(wheeling);
		    e.preventDefault();
		    var e = window.event || e; 
		    var delta = Math.max(-1, Math.min(1, (e.wheelDelta || -e.detail)));
		    
		    wheeling = setTimeout(function() {
		        wheeling = undefined;
		        if (delta === 1) {
		            moveToTop()
		        }
		        if (delta === -1) {
		            moveToBottom()
		        }
		    }, 100);

		    transition()
		}

	</script>

<script type="text/javascript" src="echart/5.0/echarts.min.js"></script>
<script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};

var option;



var data = [
  {
    name: 'High throughput data',
    children: [
      {name:'BRCA',children: [
        {name:'GSE106817 [BRCA,C=115,N=2759]',value: 1},
        {name:'GSE113486 [BRCA,C=40,N=100]',value: 1},
        {name:'GSE118782 [BRCA,C=30,N=10]',value: 1},
        {name:'GSE41922 [BRCA,C=32,N=22]',value: 1},
        {name:'GSE73002 [BRCA,C=1280,N=2686]',value: 1},
        ]},
      {name:'COAD',children: [
        {name:'GSE106817 [CRC,C=115,N=2759]',value: 1},
        {name:'GSE113486 [COAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [COAD,C=30,N=275]',value: 1},
        {name:'GSE149200 [COAD,C=5,N=5]',value: 1},
        {name:'GSE164174 [COAD,C=1423,N=1417]',value: 1},
        {name:'GSE67075 [COAD,C=18,N=30]',value: 1},
        {name:'GSE85589 [COAD,C=5,N=19]',value: 1}
        ]},
      {name:'ESCA',children: [
        {name:'GSE106817 [ESCA,C=88,N=2759]',value: 1},
        {name:'GSE112496 [ESCA,C=5,N=5]',value: 1},
        {name:'GSE112840 [ESCA,C=52,N=52]',value: 1},
        {name:'GSE113486 [ESCA,C=40,N=100]',value: 1},
        {name:'GSE122497 [ESCC,C=566,N=4965]',value: 1},
        {name:'GSE124158 [ESCA,C=30,N=275]',value: 1},
        {name:'GSE164174 [ESCA,C=1423,N=1417]',value: 1},
        {name:'GSE63108 [ESCA,C=28,N=19]',value: 1},
        ]},
      {name:'STAD',children: [
        {name:'GSE106817 [STAD,C=115,N=2759]',value: 1},
        {name:'GSE113486 [STAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [STAD,C=30,N=275]',value: 1},
        {name:'GSE126399 [STAD,C=12,N=10]',value: 1},
        {name:'GSE130654 [STAD,C=36,N=12]',value: 1},
        {name:'GSE164174 [STAD,C=1423,N=1417]',value: 1},
        {name:'GSE85589 [STAD,C=7,N=19]',value: 1}
      ]},
      {name:'HCC',children: [
        {name:'GSE106817 [HCC,C=115,N=2759]',value: 1},
        {name:'GSE113486 [HCC,C=40,N=100]',value: 1},
        {name:'GSE124158 [HCC,C=30,N=275]',value: 1},
        ]},
      {name:'LUAD',children: [
        {name:'GSE106817 [LUAD,C=115,N=2759]',value: 1},
        {name:'GSE113486 [LUAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [LUAD,C=30,N=275]',value: 1},
        {name:'GSE137140 [LUAD,C=2178,N=1746]',value: 1},
        ]},
      {name:'PAAD',children: [
        {name:'GSE106817 [PAAD,C=115,N=2759]',value: 1},
        {name:'GSE113486 [PAAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [PAAD,C=30,N=275]',value: 1},
        {name:'GSE133684 [PAAD,C=284,N=117]',value: 1},
        {name:'GSE85589 [PAAD,C=88,N=19]',value: 1},
        ]},
      {name:'BSTS',children: [
        {name:'GSE106817 [SARC,C=115,N=2759]',value: 1},
        {name:'GSE110271 [Multiple Myeloma,C=9,N=2]',value: 1},
        {name:'GSE113486 [SARC,C=40,N=100]',value: 1},
        {name:'GSE124158 [BSTS,C=892,N=275]',value: 1},
        ]},
      {name:'CHOL',children: [
        {name:'GSE113486 [CHOL,C=40,N=100]',value: 1},
        {name:'GSE144521 [CHOL,C=35,N=61]',value: 1},
        {name:'GSE85589 [CHOL,C=101,N=19]',value: 1},
        ]},
      {name:'GBM',children: [
        {name:'GSE113486 [GBM,C=40,N=100]',value: 1},
        {name:'GSE122488 [GBM,C=22,N=16]',value: 1},
        {name:'GSE124158 [GBM,C=30,N=275]',value: 1},
        {name:'GSE139031 [GBM,C=423,N=157]',value: 1},
        ]},
      {name:'PRAD',children: [
        {name:'GSE113486 [PRAD,C=40,N=100]',value: 1},
        {name:'GSE136321 [PRAD,C=24,N=24]',value: 1},
        {name:'GSE138740 [PRAD,C=146,N=89]',value: 1},
        {name:'GSE159177 [PRAD,C=278,N=187]',value: 1},
        ]},
      {name:'OV',children: [
        {name:'GSE106817 [OV,C=320,N=2759]',value: 1},
        {name:'GSE113486 [OV,C=40,N=100]',value: 1},
        {name:'GSE76449 [OV,C=24,N=4]',value: 1},
        ]},
        
      {name:'Othre',children: [
        {name:'GSE113486 [BLCA,C=392,N=100]',value: 1},
        {name:'GSE106817 [SARC,C=115,N=2759]',value: 1},
        {name:'GSE113486 [SARC,C=40,N=100]',value: 1},
        {name:'GSE110271 [Multiple Myeloma,C=9,N=2]',value: 1},
        {name:'GSE125442 [KIRC,C=10,N=10]',value: 1},
        {name:'GSE130512 [Tyriod Cancer,C=16,N=8]',value: 1},
        {name:'GSE37472 [Oral Cancer,C=50,N=26]',value: 1},
        ]},
    ]
  },
  {
    name: 'Experimental data',
    children: [
      {name: 'Homo sapiens',children: [
        {name:'lung cancer',value: 1},
        {name:'esophageal squamous cell carcinoma',value: 1},
        {name:'bladder cancer',value: 1},
        {name:'colorectal cancer',value: 1},
        {name:'multiple myeloma',value: 1},
        {name:'breast cancer',value: 1},
        {name:'liver cancer',value: 1},
        {name:'pituitary adenoma',value: 1},
        {name:'gastric cancer',value: 1},
        {name:'thyroid cancer',value: 1},
        {name:'pancreatic cancer',value: 1},
        {name:'ovarian cancer',value: 1},
        {name:'glioblastoma',value: 1},
        {name:'oral squamous cell carcinoma',value: 1},
        {name:'cervical cancer',value: 1},
        {name:'prostate cancer',value: 1},
        {name:'colon cancer',value: 1},
        {name:'glioma',value: 1},
        {name:'laryngeal cancer',value: 1},
        {name:'esophageal cancer',value: 1},
        {name:'osteosarcoma',value: 1},
        {name:'cholangiocarcinoma',value: 1},
        {name:'squamous cell carcinoma ',value: 1},
        {name:'Melanoma',value: 1},
        {name:'head and neck squamous cell carcinomas',value: 1},
        {name:'diffuse large b-cell lymphoma',value: 1},
        {name:'intracranial aneurysm ',value: 1},
        {name:'papillary thyroid cancer',value: 1},
        {name:'cholesteatoma',value: 1},
        {name:'endometrial cancer',value: 1},
        {name:'coronary artery aneurysm',value: 1},
        {name:'oral cancer',value: 1},
        {name:'acute myeloid leukemia',value: 1},
        {name:'mantle cell lymphoma',value: 1},
        {name:'chronic myeloid leukemia',value: 1},
        {name:'tongue squamous cell carcinoma',value: 1},
        {name:'urothelial cancer',value: 1},
        {name:'chronic lymphocytic leukemia',value: 1},
        {name:'nasopharyngeal carcinoma',value: 1},
        {name:'renal cancer',value: 1},
        {name:'Advanced colorectal neoplasia',value: 1},
        ]
      },
      {name: 'Mus musculus',children: [
        {name:'renal cancer',value: 1},
        {name:'liver cancer',value: 1},
        {name:'skin tumor',value: 1},
        {name:'lung cancer',value: 1},
        {name:'pancreatic cancer',value: 1},
        ]
      }
    ]
  }
];
option = {
		  title: {text: 'Data Statistics',
		 	      subtext: 'High throughput data & Experimental validation data',
				   textStyle: {
					   fontSize: 15,
					   fontFamily: 'Montserrat'
				   }},
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15}},
  series: {
    type: 'sunburst',
    data: data,
    color: ["#ff5050","#1a53ff","#ffcc00"],
    radius: [0, '90%'],
    itemStyle: {
      borderRadius: 7,
      borderWidth: 2
    },
    levels: [
      {},
      {
        r0: '15%',
        r: '35%',
        itemStyle: {
          borderWidth: 2,
          color: '#ffbf80'
        },
        label: {
          rotate: 'tangential',
          fontFamily: 'Montserrat'
        },
      },
      {
        r0: '35%',
        r: '60%',
        label: {
            fontFamily: 'Montserrat'
          },
      },
     {
        r0: '60%',
        r: '62%',
        label: {
          position: 'outside',
          padding: 3,
          silent: false,
          fontSize: 9,
          fontFamily: 'Montserrat'
        },
        itemStyle: {
          borderWidth: 3
        }
      }

    ]
  }
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>

<script type="text/javascript">
var dom = document.getElementById("container_bar");
var myChart = echarts.init(dom);
var app = {};

var option;

option = {
		  title: {text: 'Data Statistics',
	 	      subtext: 'Experimental validation data',
			   textStyle: {
				   fontSize: 15,
				   fontFamily: 'Montserrat'
			   }},
  dataset: {
    source: [
      ['score', 'amount', 'product'],
        [1,1,'Squamous cell carcinoma '],
        [1,1,'Head and neck squamous cell carcinomas'],
        [1,1,'Intracranial aneurysm '],
        [1,1,'Cholesteatoma'],
        [1,1,'Coronary artery aneurysm'],
        [1,1,'Oral cancer'],
        [1,1,'Mantle cell lymphoma'],
        [1,1,'Chronic myeloid leukemia'],
        [1,1,'Tongue squamous cell carcinoma'],
        [1,1,'Urothelial cancer'],
        [3,2,'Pituitary adenoma'],
        [3,2,'Diffuse large b-cell lymphoma'],
        [3,2,'Acute myeloid leukemia'],
        [3,2,'Chronic lymphocytic leukemia'],
        [4,2.32,'Laryngeal cancer'],
        [5,2.58,'Endometrial cancer'],
        [5,2.58,'Nasopharyngeal carcinoma'],
        [6,2.80,'Esophageal cancer'],
        [7,3,'Oral squamous cell carcinoma'],
        [8,3.16,'Osteosarcoma'],
        [9,3.32,'Skin tumor'],
        [9,3.32,'Cholangiocarcinoma'],
        [10,3.45,'Glioma'],
        [12,3.70,'Papillary thyroid cancer'],
        [12,3.70,'Advanced colorectal neoplasia'],
        [14,3.90,'Thyroid cancer'],
        [15,4,'Esophageal squamous cell carcinoma'],
        [18,4.24,'Cervical cancer'],
        [21,4.45,'Glioblastoma'],
        [28,4.85,'Ovarian cancer'],
        [31,5,'Bladder cancer'],
        [38,5.28,'Renal cancer'],
        [42,5.42,'Multiple myeloma'],
        [43,5.45,'Melanoma'],
        [52,5.72,'Gastric cancer'],
        [109,6.78,'Breast cancer'],
        [145,7.18,'Prostate cancer'],
        [156,7.29,'Lung cancer'],
        [273,8.09,'Colon cancer'],
        [585,9.19,'Colorectal cancer'],
        [757,9.56,'Pancreatic cancer'],
        [1023,10,'Liver cancer'],
    ]
  },
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15},},
  grid: { containLabel: true },
  xAxis: { name: 'count',axisLabel: {fontFamily: "Montserrat"},},
  yAxis: { type: 'category',axisLabel: {fontFamily: "Montserrat"},},
  visualMap: {
    orient: 'horizontal',
    left: 'center',
    min: 10,
    max: 100,
    textStyle:{fontFamily: 'Montserrat',fontSize:12},
    text: [' ', 'Count (Log2(X+1))'],
    // Map the score column to color
    dimension: 0,
    inRange: {
        color: ['#F5D043','#114182']
      }
  },
  series: [
    {
      type: 'bar',
        label: {
          show: true,
          precision: 2,
          position: 'right',
          valueAnimation: true,
          fontFamily: 'Montserrat'
        },
      encode: {
        // Map the "amount" column to X axis.
        x: 'amount',
        // Map the "product" column to Y axis
        y: 'product'
      }
    }
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>

</body>
</html>