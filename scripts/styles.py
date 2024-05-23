

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

    
    def table_body(self):

        body = '<tbody> \n '
        for line in self.virus_array:
            table_line = '<tr>' + ''.join([f'<td style="background-color: {self.get_color(self.pallete, i)}; width=10; heigh=5">{i}</td>' for i in line]) + '</tr>'
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
    
        
        
    

            
        
        
    


