import streamlit as st
import requests

st.set_page_config(
    page_title="حاسبة النشاط الحيوي للمركبات",
    page_icon="🧪",
    layout="centered"
)

st.title("🧪 حاسبة النشاط الحيوي للمركبات")
st.write("أدخل اسم المركب للبحث عنه في PubChem وتحليل خصائصه وحساب قاعدة Lipinski.")

st.markdown("---")

compound_name = st.text_input("اسم المركب (باللغة الإنجليزية مثل: Aspirin, Caffeine):")

if st.button("تحليل المركب", type="primary"):
    if not compound_name.strip():
        st.warning("الرجاء إدخال اسم المركب أولاً.")
    else:
        with st.spinner("جاري البحث عن المركب وجلب الخصائص من PubChem..."):
            try:
                # PubChem API Call
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name.strip()}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,Title/JSON"
                response = requests.get(url)
                
                if response.status_code == 200:
                    data = response.json()
                    props = data['PropertyTable']['Properties'][0]
                    
                    # التحويل الصريح للأرقام لتفادي أخطاء الأنواع
                    mw = float(props.get('MolecularWeight', 0))
                    logp = float(props.get('XLogP', 0))
                    hbd = int(props.get('HBondDonorCount', 0))
                    hba = int(props.get('HBondAcceptorCount', 0))
                    title = props.get('Title', compound_name)

                    st.success(f"تم العثور على المركب: {title}")
                    
                    # Display Properties
                    st.subheader("📊 الخصائص الكيميائية المستخرجة:")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("الوزن الجزئي (MW)", f"{mw} g/mol")
                        st.metric("معامل التوزيع (LogP)", logp)
                    with col2:
                        st.metric("متبرعات هيدروجين (HBD)", hbd)
                        st.metric("مستقبلات هيدروجين (HBA)", hba)

                    st.markdown("---")
                    st.subheader("⚖️ تقييم قاعدة Lipinski (Rule of Five):")

                    # Lipinski Rules Check
                    c1 = mw <= 500
                    c2 = logp <= 5
                    c3 = hbd <= 5
                    c4 = hba <= 10

                    violations = 0
                    if not c1: violations += 1
                    if not c2: violations += 1
                    if not c3: violations += 1
                    if not c4: violations += 1

                    # Results display
                    st.write(f"• الوزن الجزئي (<= 500): {'✅' if c1 else '❌'}")
                    st.write(f"• معامل LogP (<= 5): {'✅' if c2 else '❌'}")
                    st.write(f"• عدد HBD (<= 5): {'✅' if c3 else '❌'}")
                    st.write(f"• عدد HBA (<= 10): {'✅' if c4 else '❌'}")

                    if violations <= 1:
                        st.success(f"🎉 المركب يحقق شروط Lipinski! (عدد المخالفات: {violations})")
                    else:
                        st.error(f"⚠️ المركب غير مثالي كدواء طبقاً لقاعدة Lipinski (عدد المخالفات: {violations})")

                else:
                    st.error("❌ لم يتم العثور على المركب في PubChem. يرجى التأكد من كتابة الاسم باللغة الإنجليزية بشكل صحيح.")
            except Exception as e:
                st.error(f"حدث خطأ أثناء معالجة البيانات: {e}")
