mkdir -p ~/.streamlit/

echo "[general]
email = \"mwakosyajulius619@gmail.com\"
"> ~/.streamlit/credentials.toml

echo "\[server]\n\
headless = true\n\
port = $PORT\n\
enableCORS = false\n\
\n\
" > ~/.streamlit/config.toml
