<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!DOCTYPE html>
<html style="height: 100%">
<head>
<meta charset="ISO-8859-1">
<title>Insert title here</title>

</head>
<body style="height: 100%; margin: 0">
  <div id="container" style="height: 100%"></div>

  
  <script type="text/javascript" src="https://fastly.jsdelivr.net/npm/echarts@5/dist/echarts.min.js"></script>
  <!-- Uncomment this line if you want to dataTool extension
  <script type="text/javascript" src="https://fastly.jsdelivr.net/npm/echarts@5/dist/extension/dataTool.min.js"></script>
  -->
  <!-- Uncomment this line if you want to use gl extension
  <script type="text/javascript" src="https://fastly.jsdelivr.net/npm/echarts-gl@2/dist/echarts-gl.min.js"></script>
  -->
  <!-- Uncomment this line if you want to echarts-stat extension
  <script type="text/javascript" src="https://fastly.jsdelivr.net/npm/echarts-stat@latest/dist/ecStat.min.js"></script>
  -->
  <!-- Uncomment this line if you want to use map
  <script type="text/javascript" src="https://fastly.jsdelivr.net/npm/echarts@4.9.0/map/js/china.js"></script>
  <script type="text/javascript" src="https://fastly.jsdelivr.net/npm/echarts@4.9.0/map/js/world.js"></script>
  -->
  <!-- Uncomment these two lines if you want to use bmap extension
  <script type="text/javascript" src="https://api.map.baidu.com/api?v=3.0&ak=YOUR_API_KEY"></script>
  <script type="text/javascript" src="https://fastly.jsdelivr.net/npm/echarts@5/dist/extension/bmap.min.js"></script>
  -->

  <script type="text/javascript">
  
  <% 
  String cancername = "GSE142213";  
  if(request.getParameter("cancername")!=null&&!request.getParameter("cancername").equals("")){	
  	cancername = request.getParameter("cancername");
  }




  String path = request.getRealPath("/").replace("\\","/");   //get root real path 
  //String immpath  = path+"sc_cell_per/"+cancername+"_table.txt";
  String immpath  = path+"exmdb_line/qht_2_50.txt";
  String QHT = "";
  StringBuffer echart_in = new StringBuffer();
  StringBuffer echart_title = new StringBuffer();

  
  String top = "10";
  String method = "normalize";
   
  if(request.getParameter("method")!=null&&!request.getParameter("method").equals("")){	
	  method = request.getParameter("method");
// 	  System.out.println(top);
	 }
  
  Double min_value = 10.0;
  
  if(request.getParameter("top")!=null&&!request.getParameter("top").equals("")){	
	  top = request.getParameter("top");
// 	  System.out.println(top);
	 }
  
  String gseid = "BLCA_GSE145281_aPDL1";
  
  if(request.getParameter("gseid")!=null&&!request.getParameter("gseid").equals("")){	
	  gseid = request.getParameter("gseid");
// 	  System.out.println(gseid);
	 }
  
  String searchname = "GOBP_2_OXOGLUTARATE_METABOLIC_PROCESS,GOBP_3_PHOSPHOADENOSINE_5_PHOSPHOSULFATE_METABOLIC_PROCESS,GOBP_3_UTR_MEDIATED_MRNA_DESTABILIZATION,GOBP_3_UTR_MEDIATED_MRNA_STABILIZATION,GOBP_ACETYL_COA_BIOSYNTHETIC_PROCESS,GOBP_ACETYL_COA_BIOSYNTHETIC_PROCESS_FROM_PYRUVATE,GOBP_ACETYL_COA_METABOLIC_PROCESS,GOBP_ACID_SECRETION,GOBP_ACIDIC_AMINO_ACID_TRANSPORT,GOBP_ACROSOME_ASSEMBLY";

  if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
  		searchname = request.getParameter("searchname");
  		immpath  = path+"exmdb_line/"+top+"/"+gseid+"_bp.txt";
 	 }

  if(request.getParameter("searchname_cc")!=null&&!request.getParameter("searchname_cc").equals("")){	
	  	searchname = request.getParameter("searchname_cc");
// 	  	System.out.println("searchname_cc");
	  	
	  	immpath  = path+"exmdb_line/"+top+"/"+gseid+"_cc.txt";
	  }
  
  if(request.getParameter("searchname_mf")!=null&&!request.getParameter("searchname_mf").equals("")){	
	  	searchname = request.getParameter("searchname_mf");
// 	  	System.out.println("searchname_mf");
	  	
	  	immpath  = path+"exmdb_line/"+top+"/"+gseid+"_mf.txt";
	  }
  
  if(request.getParameter("searchname_hallmaker")!=null&&!request.getParameter("searchname_hallmaker").equals("")){	
	  	searchname = request.getParameter("searchname_hallmaker");
// 	  	System.out.println("searchname_hallmaker");
	  	
	  	immpath  = path+"exmdb_line/"+top+"/"+gseid+"_hallmark.txt";
	  }
  
  if(request.getParameter("searchname_pathway")!=null&&!request.getParameter("searchname_pathway").equals("")){	
	  	searchname = request.getParameter("searchname_pathway");
// 	  	System.out.println("searchname_pathway");
	  	immpath  = path+"exmdb_line/"+top+"/"+gseid+"_pathway.txt";
	  }
  
  
  
//   System.out.println(searchname);
  
  if(!searchname.equals("")){
  	for (String retval: searchname.split(",")){
  		
  		echart_title.append("'"+retval+"',");
  		min_value = dbhello.exmdb_getline_min(immpath,retval);
  		QHT = dbhello.exmdb_echart_line(immpath,retval,min_value,method);
//   	System.out.println(retval+": "+ min_value);
  		echart_in.append(QHT);
  	}	
  }

  String title_echart = echart_title.toString();
  String in_echart = echart_in.toString();

  %>
  
    var dom = document.getElementById('container');
    var myChart = echarts.init(dom, null, {
      renderer: 'canvas',
      useDirtyRect: false
    });
    var app = {};
    
    var option;

    
rawData=[["yhz","rhb","qht"],

	<%=in_echart%>

]



  // var countries = ['Australia', 'Canada', 'China', 'Cuba', 'Finland', 'France', 'Germany', 'Iceland', 'India', 'Japan', 'North Korea', 'South Korea', 'New Zealand', 'Norway', 'Poland', 'Russia', 'Turkey', 'United Kingdom', 'United States'];
  const countries = [
  <%=title_echart%>

  ];
  const datasetWithFilters = [];
  const seriesList = [];
  echarts.util.each(countries, function (rhb) {
    var datasetId = 'dataset_' + rhb;
    datasetWithFilters.push({
      id: datasetId,
      fromDatasetId: 'dataset_raw',
      transform: {
        type: 'filter',
        config: {
          and: [
            { dimension: 'qht', gte: 0 },
            { dimension: 'rhb', '=': rhb }
          ]
        }
      }
    });
    seriesList.push({
      type: 'line',
      smooth:true,
      
      datasetId: datasetId,
      showSymbol: false,
      name: rhb,
      endLabel: {
        show: true,
        fontFamily: 'Montserrat',fontSize:12,
        formatter: function (params) {
          return params.value[1] + ': ' + params.value[0];
        }
      },
      lineStyle: {width:5},
      labelLayout: {
        moveOverlap: 'shiftY'
      },
      emphasis: {
        focus: 'series',
        lineStyle: {width:10},
        itemStyle: {symbolSize: 40}
      },
      encode: {
        x: 'qht',
        y: 'Income',
        label: ['rhb', 'Income'],
        itemName: 'qht',
        tooltip: ['yhz']
      }
    });
  });
  option = {
    animationDuration: 3000,
    dataset: [
      {
        id: 'dataset_raw',
        source: rawData
      },
      ...datasetWithFilters
    ],
    title: {
      text: 'Pseudotime Pathway Analyze',
      textStyle:{fontFamily: 'Montserrat',fontSize:15},
    },
    tooltip: {
      textStyle:{fontFamily: 'Montserrat',fontSize:15},
      order: 'valueDesc',
      trigger: 'axis'
    },
    xAxis: {
      type: 'category',
      nameLocation: 'middle'
    },
    yAxis: {
      name: 'Score'
    },
    grid: {
      right: 350
    },
    series: seriesList
  };
  myChart.setOption(option);
 


    if (option && typeof option === 'object') {
      myChart.setOption(option);
    }

    window.addEventListener('resize', myChart.resize);
  </script>
</body>
</html>