import streamlit as st
import requests
import pandas as pd
import streamlit.components.v1 as components

st.set_page_config(
    page_title="منصة التحليل الكيميائي وتصميم الأدوية",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 منصة التحليل الكيميائي وتقييم الأدوية المتقدمة")
st.write("أداة شائعة لتحليل الخصائص الكيميائية، قواعد Lipinski & Veber، عرض المجسمات 3D، والمقارنة بين المركبات.")

st.markdown("---")

# إنشاء تبويبات داخل التطبيق
tab1, tab2 = st.tabs(["🔬 التحليل الفردي و 3D", "⚖️ مقارنة مركبين جنباً إلى جنب"])

# ----------------- وظيفة جلب البيانات من PubChem -----------------
def get_pubchem_data(query, search_type):
    search_by = "name" if "اسم" in search_type else "smiles"
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{search_by}/{query.strip()}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA,RotatableBondCount,Title,CID/JSON"
    res = requests.get(url)
    if res.status_code == 200:
        props = res.json()['PropertyTable']['Properties'][0]
        return {
            'cid': props.get('CID'),
            'title': props.get('Title', query),
            'mw': float(props.get('MolecularWeight', 0)),
            'logp': float(props.get('XLogP', 0)),
            'hbd': int(props.get('HBondDonorCount', 0)),
            'hba': int(props.get('HBondAcceptorCount', 0)),
            'tpsa': float(props.get('TPSA', 0)),
            'rot_bonds': int(props.get('RotatableBondCount', 0))
        }
    return None

# ================= TAB 1: التحليل الفردي =================
with tab1:
    col_in1, col_in2 = st.columns([1, 2])
    with col_in1:
        search_type = st.radio("طريقة البحث:", ("اسم المركب (Name)", "الصيغة الكيميائية (SMILES)"), key="t1_type")
    with col_in2:
        query = st.text_input("أدخل اسم المركب أو SMILES (مثال: Aspirin, Caffeine):", key="t1_query")

    if st.button("تحليل المركب", type="primary", key="btn_t1"):
        if not query.strip():
            st.warning("الرجاء إدخال البيانات أولاً.")
        else:
            with st.spinner("جاري جلب البيانات والمجسمات..."):
                data = get_pubchem_data(query, search_type)
                
                if data:
                    st.success(f"تم العثور على المركب: {data['title']} (CID: {data['cid']})")
                    
                    c_left, c_right = st.columns([1, 1])
                    
                    with c_left:
                        st.subheader("🖼️ التركيب البنائي 2D")
                        img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{data['cid']}/PNG?image_size=300x300"
                        st.image(img_url, caption=f"2D Structure: {data['title']}")

                    with c_right:
                        st.subheader("🧊 التجسيم 3D التفاعلي")
                        # Embed 3D viewer using 3Dmol.js via HTML
                        html_3d = f"""
                        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                        <div id="container-01" style="height: 300px; width: 100%; position: relative;"></div>
                        <script>
                            let viewer = $3Dmol.createViewer($("#container-01"), {{backgroundColor: "black"}});
                            $3Dmol.download("cid:{data['cid']}", viewer, {{}}, function() {{
                                viewer.setStyle({{}}, {{stick: {{}}, sphere: {{scale: 0.3}}}});
                                viewer.zoomTo();
                                viewer.render();
                            }});
                        </script>
                        """
                        components.html(html_3d, height=310)

                    st.markdown("---")
                    st.subheader("📊 الخصائص الكيميائية المستخرجة")
                    
                    m1, m2, m3, m4, m5, m6 = st.columns(6)
                    m1.metric("MW (g/mol)", data['mw'])
                    m2.metric("LogP", data['logp'])
                    m3.metric("HBD", data['hbd'])
                    m4.metric("HBA", data['hba'])
                    m5.metric("TPSA (Å²)", data['tpsa'])
                    m6.metric("Rot. Bonds", data['rot_bonds'])

                    # Lipinski & Veber
                    st.markdown("---")
                    cl1, cl2 = st.columns(2)
                    
                    with cl1:
                        st.subheader("⚖️ تقييم Lipinski")
                        c1, c2, c3, c4 = data['mw'] <= 500, data['logp'] <= 5, data['hbd'] <= 5, data['hba'] <= 10
                        v_lip = sum([not c1, not c2, not c3, not c4])
                        st.write(f"• MW <= 500: {'✅' if c1 else '❌'}")
                        st.write(f"• LogP <= 5: {'✅' if c2 else '❌'}")
                        st.write(f"• HBD <= 5: {'✅' if c3 else '❌'}")
                        st.write(f"• HBA <= 10: {'✅' if c4 else '❌'}")
                        if v_lip <= 1:
                            st.success(f"مطابق لقاعدة Lipinski (المخالفات: {v_lip})")
                        else:
                            st.error(f"غير مثالي كدواء فموي (المخالفات: {v_lip})")

                    with cl2:
                        st.subheader("📐 تقييم Veber")
                        v1, v2 = data['rot_bonds'] <= 10, data['tpsa'] <= 140
                        v_veb = sum([not v1, not v2])
                        st.write(f"• Rotatable Bonds <= 10: {'✅' if v1 else '❌'}")
                        st.write(f"• TPSA <= 140: {'✅' if v2 else '❌'}")
                        if v_veb == 0:
                            st.success("مطابق لقاعدة Veber للامتصاص المعوي الممتاز.")
                        else:
                            st.warning("قد يواجه صعوبة في الامتصاص المعوي.")

                    # CSV Export
                    st.markdown("---")
                    df_export = pd.DataFrame([data])
                    csv = df_export.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="📥 تحميل تقرير التحليل كملف CSV",
                        data=csv,
                        file_name=f"{data['title']}_analysis.csv",
                        mime='text/csv'
                    )

                else:
                    st.error("❌ لم يتم العثور على المركب في PubChem.")

# ================= TAB 2: مقارنة مركبين =================
with tab2:
    st.subheader("⚖️ قارن بين مركبين كيميائيين جنباً إلى جنب")
    
    col_m1, col_m2 = st.columns(2)
    with col_m1:
        comp1 = st.text_input("اسم المركب الأول (مثل: Aspirin):", value="Aspirin")
    with col_m2:
        comp2 = st.text_input("اسم المركب الثاني (مثل: Ibuprofen):", value="Ibuprofen")

    if st.button("بدء المقارنة", type="primary"):
        with st.spinner("جاري جلب الخصائص للمركبين..."):
            d1 = get_pubchem_data(comp1, "اسم المركب")
            d2 = get_pubchem_data(comp2, "اسم المركب")

            if d1 and d2:
                # Comparison Table
                df_comp = pd.DataFrame({
                    "الخاصية": ["اسم المركب", "الوزن الجزئي (MW)", "معامل LogP", "متبرعات (HBD)", "مستقبلات (HBA)", "المساحة القطبية (TPSA)", "الروابط القابلة للدوران"],
                    f"المركب الأول ({d1['title']})": [d1['title'], d1['mw'], d1['logp'], d1['hbd'], d1['hba'], d1['tpsa'], d1['rot_bonds']],
                    f"المركب الثاني ({d2['title']})": [d2['title'], d2['mw'], d2['logp'], d2['hbd'], d2['hba'], d2['tpsa'], d2['rot_bonds']]
                })
                
                st.table(df_comp)

                c_img1, c_img2 = st.columns(2)
                with c_img1:
                    st.image(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{d1['cid']}/PNG?image_size=250x250", caption=d1['title'])
                with c_img2:
                    st.image(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{d2['cid']}/PNG?image_size=250x250", caption=d2['title'])

            else:
                st.error("تأكد من صحة أسماء المركبات المدخلة.")
