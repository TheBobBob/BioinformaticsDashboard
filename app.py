import streamlit as st

file_path = "BioinformaticsDashboardv0.0..Rmd"

st.title("Please make sure you have RStudio and R installed before doing this!")
st.write("Click the button below to download the BioinformaticsDashboard! Currently in version 0 🧬")

with open(file_path, "rb") as file:
    st.download_button(
        label="Download File",
        data=file,
        file_name=file_path,
        mime="application/Rmd" 
    )
