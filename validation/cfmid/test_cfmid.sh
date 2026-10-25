#!/bin/bash
set -uo pipefail

MODE="${1:-msrb}"
if [ "${MODE}" != "msrb" ] && [ "${MODE}" != "cfm" ]; then
    echo "Usage: ./test_cfmid.sh [msrb|cfm]   (default: msrb)"
    exit 2
fi

IMAGE="wishartlab/cfmid:msml4.4.10_msrb1.1.15"
PLATFORM="linux/amd64"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${SCRIPT_DIR}/test_output"
INPUT_CSV="${SCRIPT_DIR}/input.csv"
OUTPUT_MSP="${WORK_DIR}/output.msp"
POS_MODEL="/trained_models_cfmid4.0/[M+H]+"
NEG_MODEL="/trained_models_cfmid4.0/[M-H]-"
POS_INPUT="${WORK_DIR}/pos_input.txt"
NEG_INPUT="${WORK_DIR}/neg_input.txt"
POS_MSP="${WORK_DIR}/pos_raw.msp"
NEG_MSP="${WORK_DIR}/neg_raw.msp"
META="${WORK_DIR}/meta.txt"
TASKS="${WORK_DIR}/tasks.tsv"
UNSUPPORTED="${WORK_DIR}/unsupported.txt"
MSRB_DIR="${WORK_DIR}/msrb"

PASS=0
FAIL=0
TOTAL=0

log() { printf '\033[1;34m[%s]\033[0m %s\n' "$(date '+%H:%M:%S')" "$*"; }
pass() { PASS=$((PASS+1)); printf '\033[1;32m  PASS\033[0m %s\n' "$*"; }
fail() { FAIL=$((FAIL+1)); printf '\033[1;31m  FAIL\033[0m %s\n' "$*"; }

log "检查 Docker 是否运行..."
if ! docker info >/dev/null 2>&1; then
    log "Docker 未运行,启动 Docker Desktop..."
    open -a Docker
    for i in $(seq 1 120); do
        if docker info >/dev/null 2>&1; then
            break
        fi
        sleep 3
    done
fi
docker info >/dev/null 2>&1 || { fail "Docker 启动失败"; exit 1; }
pass "Docker daemon 就绪"

log "拉取镜像 ${IMAGE} ..."
docker pull --platform "${PLATFORM}" "${IMAGE}" >/dev/null 2>&1 || { fail "镜像拉取失败"; exit 1; }
pass "镜像已就绪"

[ -f "${INPUT_CSV}" ] || { fail "找不到 ${INPUT_CSV}"; exit 1; }

mkdir -p "${WORK_DIR}" "${MSRB_DIR}"
rm -f "${OUTPUT_MSP}" "${POS_INPUT}" "${NEG_INPUT}" "${POS_MSP}" "${NEG_MSP}" \
    "${META}" "${TASKS}" "${UNSUPPORTED}"
rm -rf "${MSRB_DIR:?}"/*
mkdir -p "${MSRB_DIR}"

cat > "${WORK_DIR}/msrb_run.sh" <<'EOF'
#!/bin/sh
cd /data/msrb || exit 1
java -jar /cfmid/msrb-fragmenter.jar -a "$1" -ismi "$2" -o "$3" 2>&1
EOF

cat > "${WORK_DIR}/cfm_run.sh" <<'EOF'
#!/bin/sh
cfm-predict "$1" 0.001 "$2" "$3" 1 "$4" 1
EOF

log "解析 ${INPUT_CSV} (id, smiles, adduct_type, note) ..."
awk '
function csv(s, arr,    i, c, inq, n) {
    n = 1; arr[1] = ""; inq = 0
    for (i = 1; i <= length(s); i++) {
        c = substr(s, i, 1)
        if (c == "\"" && !inq) { inq = 1; continue }
        if (c == "\"" && inq) { inq = 0; continue }
        if (c == "," && !inq) { n++; arr[n] = ""; continue }
        arr[n] = arr[n] c
    }
    return n
}
function trim(x) {
    gsub(/^[ \t\r]+|[ \t\r]+$/, "", x)
    return x
}
NR>1 {
    n = csv($0, f)
    id = trim(f[1]); sm = trim(f[2]); ad = trim(f[3]); nt = trim(f[4])
    if (id == "" || sm == "") next
    print id"|"ad"|"nt
    print id"|"sm"|"ad"|"nt >> "'"${TASKS}"'"
    if (ad == "[M+H]+") print id" "sm >> "'"${POS_INPUT}"'"
    else if (ad == "[M-H]-") print id" "sm >> "'"${NEG_INPUT}"'"
    else print "id="id" adduct="ad >> "'"${UNSUPPORTED}"'"
}' "${INPUT_CSV}" > "${META}"
TOTAL=$(wc -l < "${META}" | tr -d ' ')
log "  共 ${TOTAL} 行输入, 预测引擎: ${MODE}"

if [ "${MODE}" = "cfm" ] && [ -s "${UNSUPPORTED}" ]; then
    fail "存在不支持的 adduct_type (cfm-predict 仅支持 [M+H]+ / [M-H]-):"
    while IFS= read -r line; do printf '    %s\n' "$line"; done < "${UNSUPPORTED}"
    exit 1
fi

run_group() {
    local input_file="$1"
    local model="$2"
    local out_name="$3"
    local label="$4"
    local n
    n=$(wc -l < "${input_file}" | tr -d ' ')
    log "  预测 ${label}: ${n} 个化合物"
    docker run --rm --platform "${PLATFORM}" -v "${WORK_DIR}:/data" "${IMAGE}" \
        cfm-predict "/data/${input_file##*/}" 0.001 \
        "${model}/param_output.log" "${model}/param_config.txt" \
        1 "/data/${out_name}" 1 >/dev/null 2>&1
    if [ -s "${WORK_DIR}/${out_name}" ]; then
        pass "${label} 预测完成 -> ${WORK_DIR}/${out_name}"
        return 0
    else
        fail "${label} 预测失败"
        return 1
    fi
}

annotate() {
    local raw="$1"
    local adduct="$2"
    awk -v adduct="${adduct}" '
    BEGIN {
        while ((getline line < "'"${META}"'") > 0) {
            split(line, a, "|")
            nt[a[1] "|" a[2]] = a[3]
        }
    }
    /^Name:/ {
        print
        print "Adduct: " adduct
        key = $NF "|" adduct
        if (nt[key] != "") print "Note: " nt[key]
        next
    }
    { print }
    ' "${raw}"
}

convert_msrb() {
    local src="$1" rid="$2" rad="$3" rnt="$4"
    awk -v id="${rid}" -v adduct="${rad}" -v note="${rnt}" '
    function flush() {
        if (n > 0) {
            pol = "-ve"
            if (adduct !~ /-/) pol = "+ve"
            print "Name: " pol " in-silico MS/MS by MSRB for " id
            print "Adduct: " adduct
            if (note != "") print "Note: " note
            print "ID: " id
            print "Smiles/Inchi:" sm
            print "Comment: " toupper(substr(en, 1, 1)) substr(en, 2)
            print "Num peaks: " n
            printf "%s", peaks
            print ""
        }
    }
    /^#SMILES=/ { sm = substr($0, 9); next }
    /^energy/ { flush(); en = $0; n = 0; peaks = ""; active = 1; next }
    /^$/ { active = 0; next }
    active && /^[0-9]+\.[0-9]+ [0-9]+/ {
        split($0, a, " ")
        peaks = peaks a[1] " " a[2] "\n"
        n++
    }
    END { flush() }
    ' "${src}" >> "${OUTPUT_MSP}"
}

annotate_fallback() {
    local src="$1" rid="$2" rad="$3" rnt="$4"
    awk -v id="${rid}" -v adduct="${rad}" -v note="${rnt}" '
    /^Name:/ {
        pol = "-ve"
        if (adduct !~ /-/) pol = "+ve"
        print "Name: " pol " in-silico MS/MS by CFM-ID for " id
        print "Adduct: " adduct
        if (note != "") print "Note: " note
        next
    }
    /^ID:/ { print "ID: " id; next }
    { print }
    ' "${src}" >> "${OUTPUT_MSP}"
}

> "${OUTPUT_MSP}"
OK=1

if [ "${MODE}" = "cfm" ]; then
    if [ -s "${POS_INPUT}" ]; then
        run_group "${POS_INPUT}" "${POS_MODEL}" pos_raw.msp "ESI+ [M+H]+" || OK=0
        [ "$OK" -eq 1 ] && annotate "${POS_MSP}" "[M+H]+" >> "${OUTPUT_MSP}"
    fi
    if [ -s "${NEG_INPUT}" ]; then
        run_group "${NEG_INPUT}" "${NEG_MODEL}" neg_raw.msp "ESI- [M-H]-" || OK=0
        [ "$OK" -eq 1 ] && annotate "${NEG_MSP}" "[M-H]-" >> "${OUTPUT_MSP}"
    fi
else
    n=0
    while IFS='|' read -r rid rsm rad rnt; do
        [ -n "${rid}" ] || continue
        n=$((n+1))
        log "  [${n}/${TOTAL}] ${rid} (${rad}) -> MSRB"
        docker run --rm --platform "${PLATFORM}" -v "${WORK_DIR}:/data" "${IMAGE}" \
            sh /data/msrb_run.sh "${rad}" "${rsm}" "r${n}.txt" \
            > "${MSRB_DIR}/r${n}.status" 2>&1
        msrb_out="${MSRB_DIR}/r${n}.txt_${rad}.log"
        if [ -s "${msrb_out}" ]; then
            convert_msrb "${msrb_out}" "${rid}" "${rad}" "${rnt}"
            pass "MSRB 预测完成: ${rid} -> ${msrb_out##*/}"
        else
            log "    MSRB 不支持该结构, 回退 cfm-predict"
            if printf '%s' "${rad}" | grep -q -- "-"; then
                model="${NEG_MODEL}"
            else
                model="${POS_MODEL}"
            fi
            docker run --rm --platform "${PLATFORM}" -v "${WORK_DIR}:/data" "${IMAGE}" \
                sh /data/cfm_run.sh "${rsm}" "${model}/param_output.log" \
                "${model}/param_config.txt" "/data/msrb/r${n}_fallback.msp" \
                > "${MSRB_DIR}/r${n}.fallback.status" 2>&1
            if [ -s "${MSRB_DIR}/r${n}_fallback.msp" ]; then
                annotate_fallback "${MSRB_DIR}/r${n}_fallback.msp" "${rid}" "${rad}" "${rnt}"
                pass "cfm-predict 回退完成: ${rid}"
            else
                fail "预测失败: ${rid} (${rad})"
                OK=0
            fi
        fi
    done < "${TASKS}"
fi

log "合并结果 -> ${OUTPUT_MSP}"
if [ -s "${OUTPUT_MSP}" ]; then
    records=$(grep -c "^Name:" "${OUTPUT_MSP}")
    expected=$((TOTAL * 3))
    pass "输出文件生成, 共 ${records} 条谱图记录 (预期 ${expected})"
    if [ "${records}" -eq "${expected}" ]; then
        pass "记录数校验通过: ${TOTAL} 个输入 x 3 碰撞能量"
    else
        fail "记录数不匹配: 实际 ${records}, 预期 ${expected}"
    fi
else
    fail "输出文件为空"
    OK=0
fi

log "测试完成: ${PASS} 通过, ${FAIL} 失败"
[ "$OK" -eq 1 ] && [ "${FAIL}" -eq 0 ] && exit 0 || exit 1
