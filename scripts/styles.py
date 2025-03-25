
def get_style():
    """Read CSS from external file and wrap in HTML head"""
    with open('scripts/styles.css', 'r') as f:
        css_content = f.read()
    return f"""
<head>
    <title>Sample Name</title>
    <style>
    {css_content}
    </style>
</head>
"""


class HTMLTable:
    def __init__(self, virus_2d_array, header, colorfunc, get_color):
        self.virus_array = virus_2d_array
        self.header = header
        self.pallete = colorfunc
        self.get_color = get_color

    def table_head(self):
        """Generate the table header HTML keeping the exact same output format"""
        th_elements = ''.join(f'<th>{c}</th>' for c in self.header)
        return f'<div class="tableFixHead"><table align="center">\n<thead>\n<tr>{th_elements}</tr></thead>'

    def make_td(self, x):
        if "CONTAMINATION_INFO:" in str(x):
            x, cont = x.split('CONTAMINATION_INFO:')
            return (
                '<td class="tooltip" data-tooltip="' + cont + '" '
                'style="background-color: #FFC5C5; width=10; heigh=5">' + x + '</td>'
            )
        color = self.get_color(self.pallete, x) if self.get_color else ''
        return f'<td style="background-color: {color}; width=10; heigh=5">{x}</td>'

    def table_body(self):
        """Generate table body with identical output format"""
        rows = []
        for line in self.virus_array:
            cells = ''.join(self.make_td(i) for i in line)
            rows.append(f'<tr>{cells}</tr>')
        return '<tbody>\n' + '\n'.join(rows) + '</tbody>\n</table></div>'

    def render(self):
        """Maintains the exact same method name and concatenation order"""
        return self.table_head() + self.table_body()


class HTMLDetails:
    def __init__(self, details_dict):
        self.details_dict = details_dict

    def render(self):
        """Produces identical HTML output as the original"""
        return ''.join(
            f'<details><summary>{key}</summary>{content}</details>'
            for key, content in self.details_dict.items()
        )

class VirusWarningTooltip:
    def __init__(self, ):
        self.tooltip = f'<h3 class="tooltip" style="color:#C80000">{virus_name.split("CONTAMINATION_INFO:")[0]}' \
                f'<span class="tooltip-text">{virus_name.split("CONTAMINATION_INFO:")[1]}</span></h3>'
    def __str__(self):
        return self.tooltip


class HTMLTooltip():
    def __new__(self, tooltip_text):
        return f'<span class="tooltip-text">{tooltip_text}</span>'
        

class HTMLImage:
    def __new__(self, encoded_image, type):
        self.encoded_image = encoded_image
        match type:
            case 'jpg':
                self.type = 'jpg'
            case 'svg':
                self.type = 'image/svg+xml;charset=utf-8'
        self.html_image = f'<div style="overflow: hidden; max-height:700px;"><img alt="" src="data:{self.type};base64,{self.encoded_image}" alt="Ooops! This should have been a picture" style="width: 60%; border: 2px solid #959494; min-width: 700px;"/></div> \n'
        return self.html_image



