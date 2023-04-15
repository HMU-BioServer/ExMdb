<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">






</head>



<body>
<img src="assets/images/exmdb_imm_cell/chart.png" style="position: absolute;width:28.5%;margin:1% 0">
<div id="container" style="height: 750px"></div>

<script type="text/javascript" src="echart/echarts.min.js"></script>

<script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};
var option;

<% 

String cancername = "HCC";  
if(request.getParameter("cancername")!=null&&!request.getParameter("cancername").equals("")){	
	cancername = request.getParameter("cancername");
}

String path = request.getRealPath("/").replace("\\","/");   //get root real path 
String immpath  = path+"immcell/cibersort/"+cancername+"_LM.txt";

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
	if(wh_i > 21)
	{
		wh_i = 0;
		break;
	}
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
		magicType: {type: ['stack']},
		saveAsImage: {pixelRatio: 2},
		dataView: {}
		}
  },
  tooltip: {
    trigger: 'axis',
    textStyle:{fontFamily: 'Montserrat',fontSize:15},
    position:[3,3],
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
	    type: 'value'
  },
  series: [
    {
      name: 'B_cells_naive',
      type: 'bar',
      stack: 'total',

      emphasis: {
        focus: 'series'
      },
      data: B_cells_naive
    },
    {
        name: 'B_cells_memory',
        type: 'bar',
        stack: 'total',

        emphasis: {
          focus: 'series'
        },
        data: B_cells_memory
      },
      {
          name: 'Plasma_cells',
          type: 'bar',
          stack: 'total',

          emphasis: {
            focus: 'series'
          },
          data: Plasma_cells
        },
      {
          name: 'T_cells_CD8',
          type: 'bar',
          stack: 'total',

          emphasis: {
            focus: 'series'
          },
          data: T_cells_CD8
        },
        {
            name: 'T_cells_CD4_naive',
            type: 'bar',
            stack: 'total',

            emphasis: {
              focus: 'series'
            },
            data: T_cells_CD4_naive
          },
        {
            name: 'T_cells_CD4_memory_resting',
            type: 'bar',
            stack: 'total',

            emphasis: {
              focus: 'series'
            },
            data: T_cells_CD4_memory_resting
          },
          {
              name: 'T_cells_CD4_memory_activated',
              type: 'bar',
              stack: 'total',

              emphasis: {
                focus: 'series'
              },
              data: T_cells_CD4_memory_activated
            },
            {
                name: 'T_cells_follicular_helper',
                type: 'bar',
                stack: 'total',

                emphasis: {
                  focus: 'series'
                },
                data: T_cells_follicular_helper
              },
              {
                  name: 'T_cells_regulatory_(Tregs)',
                  type: 'bar',
                  stack: 'total',

                  emphasis: {
                    focus: 'series'
                  },
                  data:T_cells_regulatory_Tregs
                },
                
                {
                    name: 'T_cells_gamma_delta',
                    type: 'bar',
                    stack: 'total',

                    emphasis: {
                      focus: 'series'
                    },
                    data:T_cells_gamma_delta
                  },
                  {
                      name: 'NK_cells_resting',
                      type: 'bar',
                      stack: 'total',

                      emphasis: {
                        focus: 'series'
                      },
                      data:NK_cells_resting
                    },
                   {
                        name: 'NK_cells_activated',
                        type: 'bar',
                        stack: 'total',

                        emphasis: {
                          focus: 'series'
                        },
                        data:NK_cells_activated
                      },
                      {
                          name: 'Monocytes',
                          type: 'bar',
                          stack: 'total',

                          emphasis: {
                            focus: 'series'
                          },
                          data:Monocytes
                        },
                        {
                            name: 'Macrophages_M0',
                            type: 'bar',
                            stack: 'total',

                            emphasis: {
                              focus: 'series'
                            },
                            data:Macrophages_M0
                          },
                          {
                              name: 'Macrophages_M1',
                              type: 'bar',
                              stack: 'total',

                              emphasis: {
                                focus: 'series'
                              },
                              data:Macrophages_M1
                            },
                            {
                                name: 'Macrophages_M2',
                                type: 'bar',
                                stack: 'total',

                                emphasis: {
                                  focus: 'series'
                                },
                                data:Macrophages_M2
                              },
                              {
                                  name: 'Dendritic_cells_resting',
                                  type: 'bar',
                                  stack: 'total',

                                  emphasis: {
                                    focus: 'series'
                                  },
                                  data:Dendritic_cells_resting
                                },
                                {
                                    name: 'Dendritic_cells_activated',
                                    type: 'bar',
                                    stack: 'total',

                                    emphasis: {
                                      focus: 'series'
                                    },
                                    data:Dendritic_cells_activated
                                  },
                                  {
                                      name: 'Mast_cells_resting',
                                      type: 'bar',
                                      stack: 'total',

                                      emphasis: {
                                        focus: 'series'
                                      },
                                      data:Mast_cells_resting
                                    },
                                    {
                                        name: 'Mast_cells_activated',
                                        type: 'bar',
                                        stack: 'total',

                                        emphasis: {
                                          focus: 'series'
                                        },
                                        data:Mast_cells_activated
                                      },
                                      {
                                          name: 'Eosinophils',
                                          type: 'bar',
                                          stack: 'total',

                                          emphasis: {
                                            focus: 'series'
                                          },
                                          data:Eosinophils
                                        },
                                        {
                                            name: 'Neutrophils',
                                            type: 'bar',
                                            stack: 'total',

                                            emphasis: {
                                              focus: 'series'
                                            },
                                            data:Neutrophils
                                          },

  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>

 

</body>
</html>