
import streamlit as st
import requests

st.set_page_config(
    page_title="حاسبة النشاط الحيوي للمركبات",
    page_icon="🧪",
    layout="centered"
)

st.title("🧪 حاسبة النشاط الحيوي وتصميم الأدوية")
st.write("ابحث باسم المركب أو صيغة SMILES لجلب الخصائص الكيميائية من PubChem وتقييم قواعد Lipinski & Veber.")

st.markdown("---")

# خيار نوع البحث
search_type = st.radio("اختر طريقة البحث:", ("اسم المركب (Name)", "الصيغة الكيميائية (SMILES)"))
query = st.text_input("أدخل المدخلات (مثال: Aspirin أو CC(=O)OC1=CC=CC=C1C(=O)O):")

if st.button("تحليل المركب الشامل", type="primary"):
    if not query.strip():
        st.warning("الرجاء إدخال اسم المركب أو صيغة SMILES أولاً.")
    else:
        with st.spinner("جاري جلب البيانات والصورة من PubChem..."):
            try:
                # تحديد الرابط بناءً على نوع البحث
                search_by = "name" if "اسم" in search_type else "smiles"
                
                # PubChem API Call
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{search_by}/{query.strip()}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA,RotatableBondCount,Title/JSON"
                response = requests.get(url)
                
                if response.status_code == 200:
                    data = response.json()
                    props = data['PropertyTable']['Properties'][0]
                    cid = props.get('CID', None)
                    
                    # تحويل الأرقام لتفادي أخطاء الأنواع
                    mw = float(props.get('MolecularWeight', 0))
                    logp = float(props.get('XLogP', 0))
                    hbd = int(props.get('HBondDonorCount', 0))
                    hba = int(props.get('HBondAcceptorCount', 0))
                    tpsa = float(props.get('TPSA', 0))
                    rot_bonds = int(props.get('RotatableBondCount', 0))
                    title = props.get('Title', query)

                    st.success(f"تم العثور على المركب: {title}")

                    # عرض الهيكل الكيميائي 2D
                    st.subheader("🖼️ التركيب البنائي للمركب (2D Structure):")
                    img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{search_by}/{query.strip()}/PNG?image_size=300x300"
                    st.image(img_url, caption=f"الشكل البنائي لـ {title}", use_column_width=False)

                    st.markdown("---")
                    
                    # Display Properties
                    st.subheader("📊 الخصائص الكيميائية المستخرجة:")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("الوزن الجزئي (MW)", f"{mw} g/mol")
                        st.metric("معامل التوزيع (LogP)", logp)
                        st.metric("المساحة القطبية (TPSA)", f"{tpsa} Å²")
                    with col2:
                        st.metric("متبرعات هيدروجين (HBD)", hbd)
                        st.metric("مستقبلات هيدروجين (HBA)", hba)
                        st.metric("روابط قابلة للدوران", rot_bonds)

                    st.markdown("---")
                    st.subheader("⚖️ 1. تقييم قاعدة Lipinski (Rule of Five):")

                    # Lipinski Rules Check
                    c1, c2, c3, c4 = mw <= 500, logp <= 5, hbd <= 5, hba <= 10
                    v_lipinski = sum([not c1, not c2, not c3, not c4])

                    st.write(f"• الوزن الجزئي (<= 500): {'✅' if c1 else '❌'}")
                    st.write(f"• معامل LogP (<= 5): {'✅' if c2 else '❌'}")
                    st.write(f"• عدد HBD (<= 5): {'✅' if c3 else '❌'}")
                    st.write(f"• عدد HBA (<= 10): {'✅' if c4 else '❌'}")

                    if v_lipinski <= 1:
                        st.success(f"🎉 مطابق لقاعدة Lipinski (عدد المخالفات: {v_lipinski})")
                    else:
                        st.error(f"⚠️ غير مثالي طبقاً لـ Lipinski (عدد المخالفات: {v_lipinski})")

                    st.markdown("---")
                    st.subheader("📐 2. تقييم قاعدة Veber (Veber's Rule):")

                    # Veber Rules Check (Rotatable Bonds <= 10 & TPSA <= 140)
                    v1, v2 = rot_bonds <= 10, tpsa <= 140
                    v_veber = sum([not v1, not v2])

                    st.write(f"• روابط قابلة للدوران (<= 10): {'✅' if v1 else '❌'}")
                    st.write(f"• المساحة السطحية القطبية TPSA (<= 140 Å²): {'✅' if v2 else '❌'}")

                    if v_veber == 0:
                        st.success("🎉 ممتاز! المركب يحقق قاعدة Veber للامتصاص المعوي الجيد.")
                    else:
                        st.warning("⚠️ المركب يخالف قاعدة Veber وقد يواجه صعوبة في الامتصاص.")

                else:
                    st.error("❌ لم يتم العثور على المركب. يرجى التأكد من كتابة الاسم أو صيغة SMILES بشكل صحيح.")
            except Exception as e:
                st.error(f"حدث خطأ أثناء معالجة البيانات: {e}")
