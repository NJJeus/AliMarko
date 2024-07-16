

set_style = """
<head>
	<title>Sample Name</title>
	 <style>
        .tableFixHead {
        overflow-y: auto; /* make the table scrollable if height is more than 200 px  */
        max-height: 500px; /* gives an initial height of 200px to the table */
        width: 95%;
        margin:auto;
        text-align: center;
                      }


        .tableFixHead thead th {
        position: sticky; /* make the table heads sticky */
        top: 0px; /* table head will be placed from the top of the table and sticks to it */
        text-align: center;
                                }


        body {
                        background-color: #f1f1f1;
                        font-family: Arial, sans-serif;
                        color: #333;
                        margin: 0;
                        padding: 0;
            text-align: center;
                }
                header {
                        background-color: #2ecc71;
                        color: #000;
                        padding: 100px;
            text-align: center;
                }
                h1 {
                        font-size: 36px;
                        margin-top: 0;
            color: #454444;
            text-align: center;
                }
        h2 {
                        font-size: 28px;
                        margin-top: 0;
            color: #454444;
            text-align: center;
        }
        th{
            z-index:1000;
        }
 hr {
  width: 800px;
  height: 4px;
  background-color: #BBB; 
  border-radius: 5px; /* smooth ends */
  border: none; /* remove default border */
  box-shadow: 0 1px 0 rgba(255, 255, 255, 0.5); /* add a subtle shadow */
}


        h3 {
        font-size: 18px;
        margin-top: 0;
        color: #454444;
        text-align: center;
                }
                main {
                        max-width: 800px;
                        margin: 0 auto;
                        padding: 10px;
            text-align: center;
                }
                p {
                        line-height: 1.5;
                        margin-bottom: 20px;
                }
                a {
                        color: #2ecc71;
                        text-decoration: none;
                }
                a:hover {
                        text-decoration: underline;
                }
                img {
                        width:95%;
                        border-radius:3px;

                }
                table {
                        width: 800px;
                        margin-bottom: 20px;
                        border-radius:4px;

            font-size: 12px;
                        box-shadow: 0px 0px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
                }
                th, td {
                        text-align: center;
                        padding: 2px;
            height: 4px;
            font-size: 17px;
                        border-bottom: 1px solid #ddd;
                }
                td {
                        background-color: #fff;

                }
                th {
                        background-color: #18aa9d;
                        color: #fff;
                        font-weight: bold;
                        border-radius:3px;

            font-size: 18px;


                }

                tr:hover {
                        background-color: #f5f5f5;
                }
                td.tooltip {
            position: relative;
        }

        td.tooltip:hover::before {
            content: attr(data-tooltip);
            position: absolute;
            background-color: #333;
            color: #fff;
            padding: 5px;
            border-radius: 6px;
            top: -30px;
            left: 50%;
            transform: translateX(-50%);
            white-space: nowrap;
            z-index: 1001;
        }
        h3.tooltip {
            position: relative;
            display: inline-block;
        }

        span.tooltip-text {
            display: none;
            position: absolute;
            background-color: #333;
            color: #fff;
            padding: 5px;
            border-radius: 6px;
            top: -30px;
            left: 100%;
            align:left;
            z-index: 1011;
            
            transform: translateX(-50%);
        }

        h3.tooltip:hover span.tooltip-text {
            display: block;
        }
        
        </style>

</head>
"""


class Table:
    def __init__(self, virus_2d_array, header, colorfunc, get_color):
        self.virus_array = virus_2d_array
        self.header = header
        self.pallete = colorfunc
        self.get_color = get_color
    def table_head(self):
        table_head = f'<div class="tableFixHead"><table align="center"> \n <thead> \n <tr>'

        
        table_head += '\n'.join([f'<th>{c}</th>' for c in self.header])
        table_head + '</tr></thead>'

        return table_head

    def make_td(self, x):
        if "CONTAMINATION_INFO:" in str(x):
            x, cont = x.split('CONTAMINATION_INFO:')
            
            return f'<td class="tooltip" data-tooltip="{cont}" style="background-color: #FFC5C5; width=10; heigh=5">{x}</td>'
        else:
            return f'<td style="background-color: {self.get_color(self.pallete, x)}; width=10; heigh=5">{x}</td>'
    
    
    def table_body(self):

        body = '<tbody> \n '
        for line in self.virus_array:
            table_line = '<tr>' + ''.join([self.make_td(i) for i in line]) + '</tr>'
            body += table_line
        body += '</tbody> \n </table></div>'
        return body
    def make_table(self):
        return self.table_head() + self.table_body()
    
class Details:
    def __init__(self, details_dict):
        self.details_dict = details_dict
    def make_details(self):
        details = ''
        for key, content in self.details_dict.items():
            details += f'<details><summary>{key}</summary>'
            details += str(content)
            details += '</details>'
        return details
    
        
        
    

            
        
        
    


