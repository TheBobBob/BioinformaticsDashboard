import streamlit as st

file_path = "BioinformaticsDashboardv0.0..Rmd"

video_path = 'BioinformaticsDashboard v0.mp4'

# Provide a description or instructions
st.title("Please make sure you have RStudio and R installed before doing this!")
st.write("Click the button below to download the BioinformaticsDashboard! Currently in version 0 🧬")

# Button to download the .Rmd file
with open(file_path, "rb") as file:
    st.download_button(
        label="Download File",
        data=file,
        file_name=file_path,
        mime="application/Rmd"  # MIME type for RMarkdown files
    )

st.write("There is a very quick demo below to show what an example output might look like! 🦠")
# Display the video
st.video(video_path)
