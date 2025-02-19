from app import app
from layout import main_layout
import callbacks  

app.layout = main_layout

if __name__ == "__main__":
    app.run_server(debug=True, host='127.0.0.1', port=8050)
