<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">






</head>



<body>
<img src="assets/images/exmdb_imm_cell/chart.png" style="position: absolute;width:28.5%;margin:1% 0">
<div id="container" style="height: 750px;"></div>

<script type="text/javascript" src="echart/echarts.min.js"></script>

<script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};
var option;

<% 

String cancername = "BRCA";  
if(request.getParameter("cancername")!=null&&!request.getParameter("cancername").equals("")){	
	cancername = request.getParameter("cancername");
}

String path = request.getRealPath("/").replace("\\","/");   //get root real path 
String immpath  = path+"immcell/xcell/Xcell_"+cancername+"_longRNAs.txt";

String cells = "";
int wh_i = 0;
BufferedReader br = new fcon(immpath).getBr();
String str = br.readLine(); // igron 1 line

String title [] = str.split("\t");



StringBuffer sb = new StringBuffer();

for(int title_i=0;title_i<title.length;title_i++){
	
	if(title_i == 0 ){
		sb.append(title[title_i]);
		sb.append("=['");
	}else if (title_i == title.length-1 ){
		sb.append(title[title_i]);
	}else{
		sb.append(title[title_i]);
		sb.append("'");
		sb.append(",");
		sb.append("'");
	}
}
sb.append("']"+"\n");

str = br.readLine();

while(str!=null && str.length()>0) {
	wh_i = wh_i + 1;
	String cellarray [] = str.split("\t");
	for(int index_i=0;index_i<cellarray.length;index_i++){
		sb.append(cellarray[index_i]);
		if(index_i == 0 ){
			sb.append("=[");
		}else{
			sb.append(",");
		}
	}
	sb.append("]"+"\n");
	str = br.readLine();
}

cells = sb.toString();
out.println(cells);
// System.out.println(cells);

%>

option = {
		  brush: {
			    toolbox: ['rect', 'polygon', 'lineX', 'lineY', 'keep', 'clear'],
			    xAxisIndex: 0
			  },
   toolbox: {
	 feature: {
		magicType: {
			type: ['stack']
			      },
			      saveAsImage: {
			          pixelRatio: 2
			        },
			      dataView: {}
			    }
			  },
  tooltip: {
    trigger: 'item',
    textStyle:{fontFamily: 'Montserrat',fontSize:15},
    axisPointer: {
      // Use axis to trigger tooltip
      type: 'shadow' // 'shadow' as default; can also be 'line' or 'shadow'
    }
  },
  legend: {
	  textStyle:{fontFamily: 'Montserrat',fontSize:15},
	  left: 'center',
	  bottom: 0,
	  type: 'scroll'
  },
  grid: {
	top:'5%',
    left: '30%',
    right: '4%',
    bottom: '10%',
    containLabel: true
  },
  xAxis: {
      type: 'category',
      data: Mixture,
      axisLabel: {
          rotate: 90,
          show:false
      }
  },
  yAxis: {
	    max:1,
	    type: 'value',
  },
  series: [
	  {name:'aDC',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: aDC},
	  {name:'Adipocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Adipocytes},
	  {name:'Astrocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Astrocytes},
	  {name:'B_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: B_cells},
	  {name:'Basophils',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Basophils},
	  {name:'CD4_memory_T_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD4_memory_T_cells},
	  {name:'CD4_naive_T_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD4_naive_T_cells},
	  {name:'CD4_T_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD4_T_cells},
	  {name:'CD4_Tcm',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD4_Tcm},
	  {name:'CD4_Tem',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD4_Tem},
	  {name:'CD8_naive_T_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD8_naive_T_cells},
	  {name:'CD8_T_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD8_T_cells},
	  {name:'CD8_Tcm',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD8_Tcm},
	  {name:'CD8_Tem',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CD8_Tem},
	  {name:'cDC',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: cDC},
	  {name:'Chondrocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Chondrocytes},
	  {name:'Class_switched_memory_B_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Class_switched_memory_B_cells},
	  {name:'CLP',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CLP},
	  {name:'CMP',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: CMP},
	  {name:'DC',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: DC},
	  {name:'Endothelial_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Endothelial_cells},
	  {name:'Eosinophils',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Eosinophils},
	  {name:'Epithelial_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Epithelial_cells},
	  {name:'Erythrocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Erythrocytes},
	  {name:'Fibroblasts',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Fibroblasts},
	  {name:'GMP',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: GMP},
	  {name:'Hepatocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Hepatocytes},
	  {name:'HSC',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: HSC},
	  {name:'iDC',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: iDC},
	  {name:'Keratinocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Keratinocytes},
	  {name:'ly_Endothelial_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: ly_Endothelial_cells},
	  {name:'Macrophages',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Macrophages},
	  {name:'Macrophages_M1',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Macrophages_M1},
	  {name:'Macrophages_M2',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Macrophages_M2},
	  {name:'Mast_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Mast_cells},
	  {name:'Megakaryocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Megakaryocytes},
	  {name:'Melanocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Melanocytes},
	  {name:'Memory_B_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Memory_B_cells},
	  {name:'MEP',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: MEP},
	  {name:'Mesangial_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Mesangial_cells},
	  {name:'Monocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Monocytes},
	  {name:'MPP',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: MPP},
	  {name:'MSC',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: MSC},
	  {name:'mv_Endothelial_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: mv_Endothelial_cells},
	  {name:'Myocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Myocytes},
	  {name:'naive_B_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: naive_B_cells},
	  {name:'Neurons',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Neurons},
	  {name:'Neutrophils',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Neutrophils},
	  {name:'NK_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: NK_cells},
	  {name:'NKT',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: NKT},
	  {name:'Osteoblast',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Osteoblast},
	  {name:'pDC',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: pDC},
	  {name:'Pericytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Pericytes},
	  {name:'Plasma_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Plasma_cells},
	  {name:'Platelets',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Platelets},
	  {name:'Preadipocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Preadipocytes},
	  {name:'pro_B_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: pro_B_cells},
	  {name:'Sebocytes',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Sebocytes},
	  {name:'Skeletal_muscle',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Skeletal_muscle},
	  {name:'Smooth_muscle',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Smooth_muscle},
	  {name:'Tgd_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Tgd_cells},
	  {name:'Th1_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Th1_cells},
	  {name:'Th2_cells',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Th2_cells},
	  {name:'Tregs',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data: Tregs},
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}


</script>

 

</body>
</html>