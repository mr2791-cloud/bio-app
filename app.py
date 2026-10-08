import streamlit as st

st.set_page_config(page_title="حاسبة النشاط الحيوي", page_icon="💊", layout="centered")

st.title("💊 حاسبة التنبؤ بالنشاط الحيوي للمركبات")
st.write("أدخل خصائص المركب الكيميائي أدناه لمعرفة ما إذا كان مثبطاً إنزيمياً واعداً أم لا:")

st.markdown("---")
mw = st.number_input("الوزن الجزيئي (Molecular Weight):", min_value=0.0, max_value=2000.0, value=300.0, step=1.0)
log_p = st.number_input("معامل الذوبان (LogP):", min_value=-10.0, max_value=10.0, value=2.5, step=0.1)

if st.button("تحليل المركب", type="primary"):
    is_mw_valid = 150 <= mw <= 400
    is_logp_valid = 1.0 <= log_p <= 4.0
    
    st.markdown("### النتيجة:")
    
    if is_mw_valid and is_logp_valid:
        st.success("المركب **[واعد وفعال]** كمثبط إنزيمي! 💊")
        st.info("السبب: الوزن الجزيئي ومعامل الذوبان في النطاق المثالي للاختراق الحيوي.")
    else:
        st.error("المركب **[غير مناسب]** أو فعاليته ضعيفة. ❌")
        st.write("**الأسباب المحتملة:**")
        if not is_mw_valid:
            st.warning(f"الوزن الجزيئي ({mw}) خارج النطاق المقبول (150 - 400).")
        if not is_logp_valid:
            st.warning(f"معامل الذوبان ({log_p}) خارج النطاق المقبول (1.0 - 4.0).")
